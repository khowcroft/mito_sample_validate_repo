import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import subprocess
import os
import pysam
import sys
import argparse
from scipy import sparse
from scipy.sparse import csc_matrix, vstack, hstack, save_npz, load_npz


# limit java?
#os.environ["_JAVA_OPTIONS"] = "-Xmx512m"

LEVEL = "0.9"
AF_FILTER = 0.9
# LEVEL = None
# AF_FILTER = None

def make_vcf_unaligned_bam(ubam_file, rcs_file, output_dir, post_filter = True, threads = 4):
    """
    Takes an UNALIGNED BAM (uBAM), aligns it to the reference using BWA, 
    sorts/indexes it, and runs Mutserve to generate a VCF.
    Inputs:
        unaligned_bam (str): Path to the input uBAM file.
        rcs_file (str): Path to the reference FASTA (e.g., rCRS.fasta).
        output_dir (str): Directory for outputs.
        post_filter (bool): If True, filters VCF for AF >= 0.9.
        threads (int): CPU threads for BWA alignment.

    Outputs:
        vcf_output (str): Path to the final VCF. Returns None on failure.
    """
    
    
    
    mt_bam = os.path.join(output_dir, "temp.bam")
    # temp = ubam_file.strip().split("/")[6]
    # print(temp)
    # mt_bam = f"/home/korb/mitochondria_sample_validation/mito_sample_validate_repo/aligned_for_ST020_and_5/{temp}_mt.bam" #remove this, I just wanted to keep the bams so I can run subsequent tests
    temp_aligned_bam = os.path.join(output_dir, "temp_aligned_sorted.bam")
    vcf_output = os.path.join(output_dir, "temp_result.vcf.gz")
    
    temp_files = [
        mt_bam, 
        mt_bam + ".bai", 
        temp_aligned_bam, 
        temp_aligned_bam + ".bai",
        vcf_output + ".tbi",
        os.path.join(output_dir, "temp_filtered.vcf.gz")
    ]
    try:
        print(f"Aligning {ubam_file} to {rcs_file} using Minimap2...")
        
        #extract fastq
        cmd_fastq = ["samtools", "fastq", "-@", str(threads), ubam_file]
        
        cmd_minimap = ["minimap2", "-ax", "lr:hq", "-t", str(threads), f"{rcs_file}.mmi", "-"]
        
        # filter out unmapped/supplementary reads BEFORE sorting
        # -u outputs uncompressed BAM (faster for piping), -F 2052 excludes unmapped/supplementary
        cmd_filter = ["samtools", "view", "-u", "-F", "2052", "-"]
        
        # sort only the successfully mapped chrM reads
        cmd_sort = ["samtools", "sort", "-@", str(threads), "-o", temp_aligned_bam, "-"]

        # Chain the processes together in memory (no intermediate hard drive writes)
        with subprocess.Popen(cmd_fastq, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as p1:
            with subprocess.Popen(cmd_minimap, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as p2:
                with subprocess.Popen(cmd_filter, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as p3:
                    subprocess.run(cmd_sort, stdin=p3.stdout, check=True)
                
        print("Indexing aligned BAM...")
        subprocess.run(["samtools", "index", temp_aligned_bam], check=True)

                
        # remove old VCF if it exists so Mutserve doesn't complain about overwriting
        if os.path.exists(vcf_output):
            os.remove(vcf_output)

        #becasue you are aligning to the rcs_file, which si only chrM
        subprocess.run([
            "samtools", "view", "-b",
            "-F", "2052",       # exclude Supplementary and unmapped
            "-o", mt_bam, 
            temp_aligned_bam 
        ], check=True)

        #index to fix warning (but works without)
        subprocess.run(["samtools", "index", mt_bam], check=True)
        coverage = get_coverage(mt_bam)


        #mutserver
        subprocess.run([
            "mutserve", "call",
            "--reference", str(rcs_file),
            "--output", vcf_output, 
            "--threads", "1",
            "--level", LEVEL, #Should adjust this number
            #"--no-freq", didnt work
            # may need to add --noFreq, because there are some that pass below --level if theyre a common variant
            mt_bam
        ], check=True)
        subprocess.run(["tabix", "-p", "vcf", vcf_output], check = True)
        
        if post_filter:
            if os.path.exists(vcf_output):
                filtered_vcf = os.path.join(output_dir, "temp_filtered.vcf.gz")
                
                # open the Mutserve output
                vcf_in = pysam.VariantFile(vcf_output)
                # open a new file for writing ('wz' writes compressed VCF)
                vcf_out = pysam.VariantFile(filtered_vcf, 'wz', header=vcf_in.header)
                
                for record in vcf_in:
                    # if a record has AF >= AF_filter, keep the record
                    keep_record = False
                    
                    for sample in record.samples.values():
                        #this will just be one sample because were only doing one sample vcf files
                        # get AF, default to 0.0 if missing
                        af = sample.get('AF', 0.0)
                        
                        # if af is a tuple or list
                        if isinstance(af, (list, tuple)):
                            af = af[0]
                        
                        if float(af) >= AF_FILTER:
                            keep_record = True
                            break
                    
                    if keep_record:
                        vcf_out.write(record)
                
                vcf_in.close()
                vcf_out.close()
                
                #overwrite the original output with the filtered one
                os.replace(filtered_vcf, vcf_output)            
            
        return vcf_output, coverage

    except subprocess.CalledProcessError as e:
        print(f"Failed to process {ubam_file}: {e}")
        return None
    
    finally:
        # Cleanup logic: Remove everything in the temp_files list
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)


def get_mito_contig_name(bam_path):
        """Inspects a BAM header to find the likely name of the mitochondrial chromosome."""
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if len(bam.references) == 1:
                return bam.references[0] # for ubam
            for ref in bam.references:
                if ref in ["chrM", "MT", "rCRS", "NC_012920.1"]:
                    return ref
        return None # If it can't find it

def make_vcf(bam_file, rcs_file, output_dir, post_filter = True):
    """
    Extracts mitochondrial reads from a BAM file, generates a VCF using Mutserve, and optionally filters it.
    
    Inputs:
        bam_file (str): Path to the input BAM file.
        rcs_file (str): Path to the reference FASTA file (e.g., rCRS.fasta).
        output_dir (str): Directory where temporary and final outputs will be saved.
        mt_chrom (str): The chromosome name for mitochondria in the BAM file (default: "chrM").
        level (str): Variant calling threshold for Mutserve (default: ".9").
        post_filter (bool): If True, filters the VCF to keep only variants with AF >= 0.9 (default: True).
        
    Outputs:
        vcf_output (str): Path to the final generated VCF file. Returns None if the process fails.
    """
    
    
    mt_chrom = get_mito_contig_name(bam_file)
    mt_bam = os.path.join(output_dir, "temp.bam")
    vcf_output = os.path.join(output_dir, "temp_result.vcf.gz")
    
    temp_files = [
        mt_bam, 
        mt_bam + ".bai", 
        vcf_output + ".tbi",
        os.path.join(output_dir, "temp_filtered.vcf.gz")
    ]
    try:
        # remove old VCF if it exists so Mutserve doesn't complain about overwriting
        if os.path.exists(vcf_output):
            os.remove(vcf_output)

        # extract Mitochondrial Reads to the temp BAM
        subprocess.run([
            "samtools", "view", "-b",
            "-F", "2052",       # exclude Supplementary and unmapped
            "-o", mt_bam, 
            str(bam_file), 
            mt_chrom
        ], check=True)

        #index to fix warning (but works without)
        subprocess.run(["samtools", "index", mt_bam], check=True)

        #mutserver
        subprocess.run([
            "mutserve", "call",
            "--reference", str(rcs_file),
            "--output", vcf_output, 
            "--threads", "1",
            "--level", LEVEL, #Should adjust this number
            #"--no-freq", didnt work
            # may need to add --noFreq, because there are some that pass below --level if theyre a common variant
            mt_bam
        ], check=True)
        subprocess.run(["tabix", "-p", "vcf", vcf_output], check = True)
        
        if post_filter:
            if os.path.exists(vcf_output):
                filtered_vcf = os.path.join(output_dir, "temp_filtered.vcf.gz")
                
                # open the Mutserve output
                vcf_in = pysam.VariantFile(vcf_output)
                # open a new file for writing ('wz' writes compressed VCF)
                vcf_out = pysam.VariantFile(filtered_vcf, 'wz', header=vcf_in.header)
                
                for record in vcf_in:
                    # if sample has AF >= af_filter, keep it
                    keep_record = False
                    
                    for sample in record.samples.values():                        
                        #this will just be one sample because were only doing one sample vcf files
                        # get AF, 0 default
                        af = sample.get('AF', 0.0)
                        
                        # if AF is tuple/list
                        if isinstance(af, (list, tuple)):
                            af = af[0]
                        
                        if float(af) >= AF_FILTER:
                            keep_record = True
                            break
                    
                    if keep_record:
                        vcf_out.write(record)
                
                vcf_in.close()
                vcf_out.close()
                
                #overwrite the original output with the filtered one
                os.replace(filtered_vcf, vcf_output)            
            

        return vcf_output

    except subprocess.CalledProcessError as e:
        print(f"Failed to process {bam_file}: {e}")
        return None
    
    finally:
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)

def get_coverage(mt_bam):
    mt_chrom = get_mito_contig_name(mt_bam)
    try:
        cov_result = subprocess.run([
                "samtools", "coverage", 
                "-r", mt_chrom, 
                "--no-header",  # Skip the header so we just get the data line
                mt_bam
            ], capture_output=True, text=True, check=True)
        
        # 1. Get the text output (e.g., "chrM  1  16569  540  100.0  55.4  ...")
        output_line = cov_result.stdout.strip()
            
        # 2. Split into a list of strings
        fields = output_line.split('\t')
            
        # 3. Extract the Mean Depth (Column 7, index 6)
        # Columns: #rname(0) start(1) end(2) numreads(3) covbases(4) coverage(5) meandepth(6) ...
        mean_depth = float(fields[6])
        return mean_depth
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools on {mt_bam}:")
        print(f"Exit code: {e.returncode}")
        print(f"Error message: {e.stderr}")
        return None 
    except (IndexError, ValueError) as e:
        print(f"Error parsing samtools output: {e}")
        return None


class MitoVariantDatabase:
    def __init__(self, n_positions=16569):
        self.cash_sum = np.array([], dtype = int) #this will make it so I don't have to sum each time (constant time?)
        # 16569 positions * 4 bases (A, C, G, T) = 66,276 possible rows
        # the empty rows in sparse matrix dont add to memory
        self.n_rows = n_positions * 4
        self.base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        
        self.db_matrix = csc_matrix((self.n_rows, 0), dtype=bool) #might change later if want to store int
        self.sample_ids = []

    def _get_row_index(self, pos, alt_base):
        """
        Calculates a fixed, 0-based row index for a specific base substitution.
        
        Inputs:
            pos (int): The 1-based genomic position from the VCF.
            alt_base (str): The alternate allele base (A, C, G, or T).
            
        Outputs:
            int: The specific row index corresponding to that mutation in the sparse matrix.
        """
        # (Pos-1 to 0-index) * 4 slots + offset for base
        return (pos - 1) * 4 + self.base_map[alt_base]

    def vcf_to_vector(self, vcf_path):
        """
        Converts a single sample's VCF into a sparse column vector representing variant presence.
        
        Inputs:
            vcf_path (str): Path to the sample's VCF file.
            
        Outputs:
            vector (scipy.sparse.csc_matrix): A sparse column vector (shape: 66276 x 1) where 1s indicate a variant.
        """
        
        vcf = pysam.VariantFile(vcf_path) #returns an iterator ignoring the headers and unziping if needed
        active_rows = []
        for record in vcf:
            #skip if not a snp (change this later?)
            if len(record.ref) != 1: 
                continue 
            
            for alt in record.alts:
                if len(alt) != 1: #skips insertion
                    continue
                
                
                try: #get the row index
                    idx = self._get_row_index(record.pos, alt)
                    if 0 <= idx < self.n_rows:
                        active_rows.append(idx)
                except KeyError:
                    # Handles cases like 'N' or non-standard bases
                    continue
        
        # data and cols doesn't really matter, but still
        # might need to change these to ints later
        # just keep as ints for now, (this is like a MB then?)
        data = np.ones(len(active_rows), dtype=int)
        rows = np.array(active_rows)
        cols = np.zeros(len(active_rows), dtype=int)
        
        vector = sparse.coo_matrix(
            (data, (rows, cols)), 
            shape=(self.n_rows, 1)
        ).tocsc()
        
        return vector

    def add_sample(self, vcf_path, sample_id, personal_id, coverage=None):
        """
        Adds a new sample to the existing sparse matrix database.
        
        Inputs:
            vcf_path (str): Path to the VCF file to be added.
            sample_id (str): Unique identifier for the acquisition/sample.
            personal_id (str): Identifier for the biological source.
            coverage (float): Mean sequencing depth of the sample.
            
        Outputs:
            None. Modifies self.db_matrix, self.sample_ids, and self.cash_sum in place.
        """
        if ([sample_id, str(personal_id), coverage] in self.sample_ids):
            print("sample already in database, not actually adding")
        else:
            new_vector = self.vcf_to_vector(vcf_path)
            # add the new column to existing DB
            self.db_matrix = hstack([self.db_matrix, new_vector])
            self.sample_ids.append((sample_id, str(personal_id), coverage))
            
            self.cash_sum = np.append(self.cash_sum, new_vector.sum())
        
    def remove_sample(self, target_sample_id):
        """
        Removes a sample from the database using its sample_id.
        
        Inputs:
            target_sample_id (str): The unique sample_id (acq_id) to remove.
            
        Outputs:
            bool: True if removed successfully, False if the sample was not found.
        """
        # get sample index
        target_idx = -1
        for i, (s_id, b_id, cov) in enumerate(self.sample_ids):
            if s_id == target_sample_id:
                target_idx = i
                break
                
        if target_idx == -1:
            print(f"Sample '{target_sample_id}' not found in the database.")
            return False
            
        # remove from json
        removed_metadata = self.sample_ids.pop(target_idx)
        
        # remove from matrix
        mask = np.ones(self.db_matrix.shape[1], dtype=bool)
        mask[target_idx] = False
        self.db_matrix = self.db_matrix[:, mask]
        
        # remove from sm
        self.cash_sum = np.delete(self.cash_sum, target_idx)
        
        print(f"Successfully removed sample: {removed_metadata[0]} (PID: {removed_metadata[1]})")
        return True

    def compare_sample(self, vcf_path, threshold=0.9):
        """
        Compares a new sample against all samples in the database using dot-product intersection.
        Formula: sum(dot) / min(sum_db, sum_sample)
        might add ISB0
        Note: if there is not a variant it will just say it's different from every sample

        Inputs:
            vcf_path (str): Path to the VCF file to compare.
            threshold (float): The similarity ratio required to consider a sample a match (default: 0.9).
            
        Outputs:
            numpy.ndarray (bool): A boolean array indicating which database samples meet the match threshold.
        """
        if self.db_matrix.shape[1] == 0:
            return np.array([])
        
        if isinstance(vcf_path, str):
        # Handle string or path string
            # print(f"Variable is a string")
            sample_vec = self.vcf_to_vector(vcf_path)
            # print("sample_vec has shape: ", sample_vec.shape)
        elif isinstance(vcf_path, csc_matrix):
            # print(f"Variable is csc_matrix wiht shape {vcf_path.shape}")
            sample_vec = vcf_path
        else:
            print("Unknown type")

        
        
        #1: dot product of sample with each other sample and summed
        intersection = self.db_matrix.T @ sample_vec
        intersection = intersection.toarray().flatten()

        #2: denominator is the min of each sample comparison
        db_sums = self.cash_sum
        sample_sum = sample_vec.sum()
        
        denominators = np.minimum(db_sums, sample_sum)
        
        #3: this is a more robust way of just doing the division in case there is a sample that matches perfectly
        #or doesn't have the reads needed it just will put a 0, meaning different from every files
        with np.errstate(divide='ignore', invalid='ignore'):
            scores = intersection / denominators
            scores = np.nan_to_num(scores)

        return scores > threshold
    
    
    
    def visualize_sample(self, vcf_path, outfile, coverage = 100.0, threshold=0.8, compare = None, target_sample = None): #set to 100 so it works so the default isn't low coverage
        """
        #get the database where the columns are the ones that match
        
        #make a plot like igv, where the the length is 16k bases, and each snp in sample_vec is shown.
        #I want each of the matching samples to be shown above with the matching snps highlighted green and non matching ones red
        """
        
        """
        Generates an IGV-style plot comparing a target sample's variants against matching database samples.
        
        Inputs:
            vcf_path (str): Path to the target VCF file or array if you already have it
            coverage (float): Coverage of the target sample; influences plotting color (default: 100.0).
            threshold (float): Minimum similarity threshold to retrieve matches for plotting (default: 0.8).
            compare (str): Specific PID to force a comparison with, overriding the threshold (default: None).
            
        Outputs:
            None. saves a matplotlib plot.
        """
        
        if isinstance(vcf_path, str):
        # Handle string or path string
            # print(f"Variable is a string")
            sample_vec = self.vcf_to_vector(vcf_path)
            # print("sample_vec has shape: ", sample_vec.shape)
        elif isinstance(vcf_path, csc_matrix):
            # print(f"Variable is csc_matrix wiht shape {vcf_path.shape}")
            sample_vec = vcf_path
        else:
            print("Unknown type")
            print(vcf_path.type())
        

        # # 1. Parse the input sample
        # sample_vec = self.vcf_to_vector(vcf_path)
        # Get indices of variants (rows in the sparse matrix)
        target_indices = sample_vec.indices 
        # Convert row indices back to genome position (1-16569)
        target_pos = (target_indices // 4) + 1
        
        # 2. Find matches
        #print(compare)
        if compare:
            match_idxs = []
            for i, sample in enumerate(self.sample_ids):
                if sample[1] == compare:
                    match_idxs.append(i)
        else: 
            matches_mask = self.compare_sample(vcf_path, threshold)
            match_idxs = np.where(matches_mask)[0]

        # print(len(match_idxs))
        # print(match_idxs)
        if len(match_idxs) == 0:
            print(f"No samples found with similarity > {threshold}")
            return

        print(f"Visualizing {len(match_idxs)} matching samples...")

        # 3. Setup Plot
        # Dynamic height based on number of matches
        fig, ax = plt.subplots(figsize=(14, 2 + 0.5 * len(match_idxs)))
        
        y_ticks = [0]
        if target_sample:
            y_labels = [target_sample]
        else: y_labels = ["Target Input"]

        # 4. Plot Target Sample (Reference for this view) at y=0
        print(coverage)
        if coverage > 10:
            ax.scatter(target_pos, np.zeros_like(target_pos), 
                c='blue', marker='|', s=100, label='Target Variants', zorder=10)
        else:
            ax.scatter(target_pos, np.zeros_like(target_pos), 
                c='orange', marker='|', s=100, label='Target Variants', zorder=10)

        # 5. Plot Matching Samples stacked above
        for i, db_idx in enumerate(match_idxs):
            y_pos = i + 1
            
            # Retrieve Metadata
            s_id, b_id, coverage_comp = self.sample_ids[db_idx]
            
            # Retrieve Sparse Vector for this DB sample
            db_col = self.db_matrix[:, db_idx]
            db_indices = db_col.indices
            
            # --- logic for coloring ---
            # Intersection (Shared variants)
            shared_indices = np.intersect1d(db_indices, target_indices)
            # Difference (Variants in DB but not in Target)
            unique_db_indices = np.setdiff1d(db_indices, target_indices)
            
            # Convert to positions
            shared_pos = (shared_indices // 4) + 1
            unique_pos = (unique_db_indices // 4) + 1
            
            # Plot Green (Matches)
            if coverage_comp > 10:
                if len(shared_pos) > 0:
                    ax.scatter(shared_pos, np.full_like(shared_pos, y_pos), 
                            c='green', marker='o', s=30, alpha=0.8, edgecolors='none')
                
                # Plot Red (Mismatches)
                if len(unique_pos) > 0:
                    ax.scatter(unique_pos, np.full_like(unique_pos, y_pos), 
                            c='red', marker='x', s=30, alpha=0.8)
            else:
                if len(shared_pos) > 0:
                    ax.scatter(shared_pos, np.full_like(shared_pos, y_pos), 
                            c='orange', marker='o', s=30, alpha=0.8, edgecolors='none')
                
                # Plot Red (Mismatches)
                if len(unique_pos) > 0:
                    ax.scatter(unique_pos, np.full_like(unique_pos, y_pos), 
                            c='orange', marker='x', s=30, alpha=0.8)

            y_ticks.append(y_pos)
            y_labels.append(f"{s_id} ({b_id})")

        # 6. Formatting
        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_labels)
        ax.set_xlabel("Mitochondrial Position (bp)")
        ax.set_title(f"Variants Comparison")
        ax.set_xlim(0, 16569)
        ax.set_ylim(-1, len(match_idxs) + 1)
        
        # Add grid lines for easier reading of positions
        ax.grid(True, axis='x', linestyle=':', alpha=0.5)
        
        # Custom Legend
        legend_elements = [
            Line2D([0], [0], color='blue', marker='|', linestyle='None', markersize=10, label='Target Variants'),
            Line2D([0], [0], color='orange', marker='|', linestyle='None', markersize=10, label='Target Variants (low coverage sample)'),
            Line2D([0], [0], color='green', marker='o', linestyle='None', markersize=8, label='Shared (Match)'),
            Line2D([0], [0], color='red', marker='x', linestyle='None', markersize=8, label='Mismatch (DB Only)'),
            Line2D([0], [0], color='orange', marker='o', linestyle='None', markersize=8, label='Shared (Match), low coverage sample'),
            Line2D([0], [0], color='orange', marker='x', linestyle='None', markersize=8, label='Mismatch (DB Only), low coverage sample')

        ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))

        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()
        
    def compare_for_hist(self, vcf_path):
        """
        Calculates the raw similarity scores between a target sample and the entire database.
        
        Inputs:
            vcf_path (str): Path to the VCF file to evaluate.
            
        Outputs:
            numpy.ndarray (float): An array of continuous similarity scores for each database sample.
        """
        if self.db_matrix.shape[1] == 0:
            return np.array([])

        sample_vec = self.vcf_to_vector(vcf_path)
        
        #1: dot product of sample with each other sample and summed
        intersection = self.db_matrix.T @ sample_vec
        intersection = intersection.toarray().flatten()

        #2: denominator is the min of each sample comparison
        db_sums = np.array(self.db_matrix.sum(axis=0)).flatten()
        sample_sum = sample_vec.sum()
        
        denominators = np.minimum(db_sums, sample_sum)
        
        #3: this is a more robust way of just doing the division in case there is a sample that matches perfectly
        #or doesn't have the reads needed it just will put a 0, meaning different from every files
        with np.errstate(divide='ignore', invalid='ignore'):
            scores = intersection / denominators
            scores = np.nan_to_num(scores)

        return scores
    
    def make_output(self, matches, mito_db, acq_id, personal_id, output_dir, file_name = "matches_log.txt"):
        """
        Appends matching samples to a text file.
        
        Inputs:
            matches (numpy.ndarray): Boolean array indicating matched samples.
            mito_db (MitoVariantDatabase): The database instance containing sample metadata.
            acq_id (str): Acquisition ID of the evaluated sample.
            personal_id (str): PID of the evaluated sample.
            output_dir (str): Directory where the log file will be stored.
            file_name (str): Name of the output log file (default: "matches_log.txt").
            
        Outputs:
            None. Writes out to the filesystem.
        """
        if np.any(matches):
            # 1. Get the row numbers
            match_indices = np.where(matches)[0]
            
            # 2. Convert row numbers to Sample Names (e.g., [0, 5] -> ["Sample_A", "Sample_F"])
            matching_names = [mito_db.sample_ids[i] for i in match_indices]
            
            # 3. Create the message string
            log_message = f"Sample {acq_id}, {personal_id} matches existing samples: {matching_names}\n"
                        
            # 5. Append to File
            # mode='a' means "append" (add to the end without deleting previous lines)
            log_path = os.path.join(output_dir, file_name)
            with open(log_path, "a") as f:
                f.write(log_message)
                
                

        

    def create_database_from_csv(self, rcs_file, output_dir, csv_file, save_file):
        """
        Iterates over a CSV to process BAM files, extract variants, and build a new database from scratch.
        
        Inputs:
            rcs_file (str): Path to the reference FASTA file.
            output_dir (str): Directory for temporary processing files.
            csv_file (str): Path to the input CSV mapping acquisitions to IDs.
            save_file (str): File prefix where the resulting .npz and .json files will be saved.
            
        Outputs:
            None. Saves the new database state to disk.
        """
        df = pd.read_csv(csv_file, index_col = 0)
        filtered_df = df[df["acquisitions"].notna()]

        for index, row in filtered_df.iterrows():
            acq_folder = row['acquisitions']
            acq_id = row['acq_id']
            extension = row['extension']
            sample_bam = acq_folder + f"/{acq_id}{extension}"
            
            personal_id = str(row.name)
            
            if is_unaligned_bam(sample_bam):
                vcf_file, coverage = make_vcf_unaligned_bam(sample_bam, rcs_file, output_dir)
            else:
                vcf_file = make_vcf(sample_bam, rcs_file, output_dir)
                coverage = get_coverage(sample_bam)
            
            
            #check matches
            #matches = mito_db.compare_sample(vcf_file, threshold=THRESHOLD) 
            
            #make_output(matches, mito_db, acq_id, personal_id, output_dir)

            # add to Database
            self.add_sample(vcf_file, sample_id=acq_id, personal_id = personal_id, coverage = coverage)

        #save Final Database
        self.save(save_file)
        
    def make_igv_output_putative(self, putative_set, threshold, output_dir):
        """helper to make the output of every putative error into a file for checkall

        Args:
            putative_set (set of tuples): a set of the indices of the rows to look at and whether it's putative due to low coverage
            threshold (float): the threshold at which to consider files the same PID
            output_dir (str): the output directory
        """
        putative_dir = output_dir+"/putative_visualization"
        coverage_dir = putative_dir+"/low_coverage"
        os.makedirs(putative_dir, exist_ok=True)
        os.makedirs(coverage_dir, exist_ok=True)
        # print(putative_set)
        # assert(False)
        for index, low_coverage in putative_set:
            array_1 = self.db_matrix[:, index]
            coverage = self.sample_ids[index][2]
            
            #make a visualization of sample 1 against its supposed sample
            coverage_str = ""
            if low_coverage:
                coverage_str = "low_coverage/"
            first_file = f"{putative_dir}/{coverage_str}compare_{self.sample_ids[index][0]}_to_sample_{self.sample_ids[index][1]}"
            self.visualize_sample(vcf_path = array_1, outfile= first_file, coverage=coverage, threshold = threshold, compare = self.sample_ids[index][1], target_sample=self.sample_ids[index][0])
            
            #make vis against everything that matches
            second_file = f"{putative_dir}/{coverage_str}compare_{self.sample_ids[index][0]}_to_matches"
            self.visualize_sample(vcf_path =array_1, outfile= second_file, coverage=coverage, threshold = threshold, target_sample=self.sample_ids[index][0])
            
    def make_igv_output_putative_compare(self, threshold, output_dir, filter_sample = None):
        
        putative_dir = output_dir+"/putative_visualization"+f"/{filter_sample}"
        os.makedirs(putative_dir, exist_ok=True)

        # Find the index
        index = None
        for i, sublist in enumerate(self.sample_ids):
            if sublist[0] == filter_sample:
                index = i
                break
        
        array_1 = self.db_matrix[:, index]
        coverage = self.sample_ids[index][2]
            
        #make a visualization of sample 1 against its supposed sample
            
        first_file = f"{putative_dir}/compare_{self.sample_ids[index][0]}_to_sample_{self.sample_ids[index][1]}"
        self.visualize_sample(vcf_path = array_1, outfile= first_file, coverage=coverage, threshold = threshold, compare = self.sample_ids[index][1], target_sample=self.sample_ids[index][0])
            
        #make vis against everything that matches
        second_file = f"{putative_dir}/compare_{self.sample_ids[index][0]}_to_matches"
        # print(second_file)
        self.visualize_sample(vcf_path =array_1, outfile= second_file, coverage=coverage, threshold = threshold, target_sample=self.sample_ids[index][0])
            


    def check_all_database(self, output_dir, histogram_save, log_file, threshold, visualize_putative = False):
        """
        Calculates pairwise similarity for all samples currently in the database 
        and plots a histogram of the scores. Excludes self-comparisons.
        

        Args:
            output_dir (str): the destination of the histogram, log file, and visualizations
            histogram_save (str): the name of the histogram file
            log_file (str): the name of the log_file
            threshold (float): the threshold at which to consider two files the same or different (default is .8)
            visualize_putative (bool, optional): when this is true, you get a folder of visualizations for everything that is potentially incorrect based on the given labels. Defaults to False.
            filter_putative (str, optional): this is useful for when you add to a database and want to look at the potentail errors from exactly one same rather than the whole database. Defaults to None.
        """
        n_samples = self.db_matrix.shape[1]
        if n_samples < 2:
            print("Database needs at least 2 samples for a pairwise comparison.")
            return

        print(f"Calculating pairwise similarities for {n_samples} samples...")

        # dot product of all columns against all columns
        # gives NxN where each column is the comparisons to a sample (i guess each row is too)
        intersection_matrix = (self.db_matrix.T @ self.db_matrix).toarray()

        # min(sum_i, sum_j) for all possible pairs
        # takes the sums [1,2,3], makes [1,2,3] and [1],[2],[3] and then
        #broadcasts into the 3x3 to take the minimum of each comparison
        sums = self.cash_sum
        denominators = np.minimum(sums[:, None], sums[None, :])

        # scores
        with np.errstate(divide='ignore', invalid='ignore'):
            scores_matrix = intersection_matrix / denominators
            scores_matrix = np.nan_to_num(scores_matrix)

        # upper triangle to avoid self-duplicates
        upper_tri_idx = np.triu_indices(n_samples, k=1)
        all_scores = scores_matrix[upper_tri_idx]

        # log pairs
        log_path = os.path.join(output_dir, log_file)
        print(f"Logging to {log_path}...")
        
        with open(log_path, "w") as f:
            visaully_check_set = set()
            for i, j in zip(upper_tri_idx[0], upper_tri_idx[1]):
                coverage_warning = ""
                cov_end = ""
                low_cov_i_or_j = False
                if self.sample_ids[i][2] < 10:
                    coverage_warning += f"Warning: coverage on sample {self.sample_ids[i][0]} = {self.sample_ids[i][2]} < 10\n"
                    cov_end = "\n"
                    low_cov_i_or_j = True
                if  self.sample_ids[j][2] < 10:
                    coverage_warning += f"Warning: coverage on sample {self.sample_ids[j][0]} = {self.sample_ids[j][2]} < 10\n"
                    cov_end = "\n"
                    low_cov_i_or_j = True
                if scores_matrix[i, j] >= threshold:
                    if self.sample_ids[i][1] != self.sample_ids[j][1]:
                        # sample_ids elements are tuples like (sample_id, personal_id, coverage)
                        s1 = f"{self.sample_ids[i][0]} ({self.sample_ids[i][1]})"
                        s2 = f"{self.sample_ids[j][0]} ({self.sample_ids[j][1]})"
                        f.write(f"{cov_end}{coverage_warning}Match above threshold for different PIDs ({scores_matrix[i, j]:.2f}): {s1} <--> {s2}{cov_end}\n")
                        
                        #making the igv esque output for each failing sample:
                        if visualize_putative:
                            visaully_check_set.add((i, low_cov_i_or_j))
                            visaully_check_set.add((j, low_cov_i_or_j))
                    elif scores_matrix[i, j] < 1:
                        s1 = f"{self.sample_ids[i][0]} ({self.sample_ids[i][1]})"
                        s2 = f"{self.sample_ids[j][0]} ({self.sample_ids[j][1]})"
                        f.write(f"{cov_end}{coverage_warning}Match above threshold but below 1 same PID ({scores_matrix[i, j]:.2f}): {s1} <--> {s2}{cov_end}\n")
                elif self.sample_ids[i][1] == self.sample_ids[j][1]:
                    s1 = f"{self.sample_ids[i][0]} ({self.sample_ids[i][1]})"
                    s2 = f"{self.sample_ids[j][0]} ({self.sample_ids[j][1]})"
                    f.write(f"{cov_end}{coverage_warning}Match below threshold from same PID ({scores_matrix[i, j]:.2f}): {s1} <--> {s2}{cov_end}\n")
                    if visualize_putative:
                        visaully_check_set.add((i, low_cov_i_or_j))
                        visaully_check_set.add((j, low_cov_i_or_j))
                            
        if visualize_putative:
            self.make_igv_output_putative(visaully_check_set, threshold, output_dir)

        #plot
        # Extract personal IDs and coverages into arrays for quick comparison
        pids = np.array([s[1] for s in self.sample_ids])
        covs = np.array([s[2] for s in self.sample_ids])
        
        # Get the PIDs and coverages for the pairs we are plotting
        pids_i = pids[upper_tri_idx[0]]
        pids_j = pids[upper_tri_idx[1]]
        covs_i = covs[upper_tri_idx[0]]
        covs_j = covs[upper_tri_idx[1]]
        
        # Create boolean masks for same vs different PIDs
        same_id_mask = (pids_i == pids_j)
        diff_id_mask = ~same_id_mask

        # Create boolean masks for coverage (True if EITHER sample in the pair is < 10)
        low_cov_mask = (covs_i < 10) | (covs_j < 10)
        high_cov_mask = ~low_cov_mask
        
        # Filter scores into the four new categories
        scores_same_high = all_scores[same_id_mask & high_cov_mask]
        scores_diff_high = all_scores[diff_id_mask & high_cov_mask]
        scores_same_low = all_scores[same_id_mask & low_cov_mask]
        scores_diff_low = all_scores[diff_id_mask & low_cov_mask]
            
        plt.figure(figsize=(10, 6))
            
        # Plot stacked using the 4 categories
        plt.hist([scores_same_high, scores_diff_high, scores_same_low, scores_diff_low], 
                bins=50, 
                stacked=True, 
                color=['#2ca02c', '#d62728', 'yellowgreen', 'tomato'], 
                label=[
                    'Same PID (Good Coverage)', 
                    'Different PID (Good Coverage)',
                    'Same PID (Low Coverage)', 
                    'Different PID (Low Coverage)'
                ],
                edgecolor='black',
                alpha=0.7)
        
        plt.axvline(threshold, color='black', linestyle='dashed', linewidth=1.5, label=f'Threshold ({threshold})')
            
        plt.title(f"Pairwise Score Distribution (N={len(all_scores)} comparisons)")
        plt.xlabel("Similarity Score")
        plt.ylabel("Frequency")
        plt.legend()
            
        save_path = os.path.join(output_dir, histogram_save)
        plt.savefig(save_path)
        print(f"Histogram saved to {save_path}")
        plt.show()
    
    def save(self, output_prefix):
        """
        Serializes the sparse matrix and sample metadata to disk.
        
        Inputs:
            output_prefix (str): Prefix path for the .npz (matrix) and .json (metadata) files.
            
        Outputs:
            None.
        """
        if self.db_matrix is not None:
            save_npz(f"{output_prefix}.npz", self.db_matrix)
            
            #save sample_ids
            import json
            with open(f"{output_prefix}_meta.json", "w") as f:
                json.dump({
                    "sample_ids": self.sample_ids
                }, f)
                
    def load(self, input_prefix):
        """
        Deserializes and loads the sparse matrix and metadata from disk into the current instance.
        
        Inputs:
            input_prefix (str): Prefix path of the stored .npz and .json files.
            
        Outputs:
            None. Replaces internal state in place.
        """
        import json
        
        # load the Matrix
        # this reconstructs the CSC matrix object from the binary file
        self.db_matrix = load_npz(f"{input_prefix}.npz")
        
        self.cash_sum = np.array(self.db_matrix.sum(axis=0)).flatten()
        
        # load the Metadata
        with open(f"{input_prefix}_meta.json", "r") as f:
            meta = json.load(f)
            self.sample_ids = meta["sample_ids"]
            
        print(f"Database loaded: {self.db_matrix.shape[1]} samples.")



def remove_files(outdir):
    temp_file_path = os.path.join(outdir, "temp_result.txt")
    if os.path.exists(temp_file_path):
        os.remove(temp_file_path)
    temp_file_path = os.path.join(outdir, "temp_result.vcf.gz")
    if os.path.exists(temp_file_path):
        os.remove(temp_file_path)
        
    index_path = os.path.join(outdir, "temp_result.vcf.gz.tbi")
    if os.path.exists(index_path):
        os.remove(index_path)
        
def is_unaligned_bam(bam_path):
    """
    Checks if a BAM file is unaligned (uBAM) or aligned.
    Returns True if unaligned, False if aligned.
    """
    try:
        # (check_sq=False allows reading uBAMs)
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
            # 1. if there are no reference sequences defined, it's definitely a uBAM.
            if len(bam.references) == 0:
                return True
            
            # 2. if references exist, check the first few reads.
            # IF FIRST 25 reads unmapped assume ubam. 
            unmapped_count = 0
            check_limit = 25
            
            for i, read in enumerate(bam):
                if i >= check_limit:
                    break
                if read.is_unmapped:
                    unmapped_count += 1
                    
            if unmapped_count == check_limit:
                return True
                
            return False
            
    except Exception as e:
        print(f"Error reading BAM file {bam_path}: {e}")
        # attempt alignment if something goes weird
        return True
        
        
def get_or_create_mmi(rcs_file):
    """
    Checks if a Minimap2 index (.mmi) exists for the given reference FASTA.
    If not, it generates one.
    
    Inputs:
        rcs_file (str): Path to the reference FASTA file.
    Outputs:
        mmi_file (str): Path to the Minimap2 index file.
    """
    mmi_file = f"{rcs_file}.mmi"
    
    if os.path.exists(mmi_file):
        print(f"Found existing Minimap2 index: {mmi_file}")
    else:
        print(f"Minimap2 index not found. Building {mmi_file}...")
        try:
            subprocess.run(["minimap2", "-d", mmi_file, rcs_file], check=True)
            print("Index built successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to build Minimap2 index: {e}")
            
def main():
    parser = argparse.ArgumentParser(description="Mitochondrial Variant Database Manager")
    
    # Main actions (mutually exclusive)
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument('-add', action='store_true', help='Add a sample to an existing database')
    action_group.add_argument('-remove', action='store_true', help='Remove a sample from an existing database')
    action_group.add_argument('-create', action='store_true', help='Create a database from a CSV file')
    action_group.add_argument('-compare', action='store_true', help='Compare a sample with the database')
    action_group.add_argument('-checkall', action='store_true', help='Compare all samples within the entire database')

    # Other important
    parser.add_argument('--outdir', type=str, help='Output directory for logs and temp files')
    parser.add_argument('--db', type=str, required=True, help='Path prefix for the database (e.g., /path/to/my_db) (as input or output)')
    parser.add_argument('--rcs', type=str, help='Path to reference FASTA for chrM only (Required for -add, -create, -compare)')
    parser.add_argument('--threshold', type=float, default=0.8, help='Similarity threshold (default: 0.8)')
    
    # main helpers
    parser.add_argument('--csv', type=str, help='Path to the CSV file (Required for -create)')
    parser.add_argument('--bam', type=str, help='Path to the BAM file (Required for -add and -compare)')
    parser.add_argument('--acq', type=str, help='Acquisition ID (Required for -add and -remove) This is the individual sample ID')
    parser.add_argument('--pid', type=str, help='personal ID (Required for -add) This is the group ID, which can contain multiple acquisition IDs')
    parser.add_argument('--hist', type=str, default='score_histogram.png', help='Histogram output filename (Used with -checkall)')
    parser.add_argument('--log_file', type=str, default='pairwise_matches_log.txt', help='Name of log file (Used with -checkall)')
    parser.add_argument('--visualize', type=str, default='visualization.png', help='Visualize file matches after comparing (Optional for -compare)')
    parser.add_argument('--compare_sample', type=str, default = "", help='A specific PID to compare a sample to (not just the ones it matches) (Optional for -compare)')
    parser.add_argument('--vis_putative', action='store_true', help='when checking against the entire database, it will show show a visualization of all the potentially erroneous samples for visual inspection')
    parser.add_argument('--filter_sample', type=str, default = "", help='when doing compare this allows you to get any putative errors related to a particular sample (sample must be in database already, otherwise put the bam in --bam)')



    #temp to test levels
    # global LEVEL, AF_FILTER
    # parser.add_argument('--af_filter', type=float, default=0.9, help='test')
    # parser.add_argument('--level', type=str, default=".9", help='test')
    
    # args = parser.parse_args()
    
    # # assign the parsed values to the globals
    # LEVEL = args.level
    # AF_FILTER = args.af_filter

    # print(LEVEL)
    # print(AF_FILTER)
    args = parser.parse_args()
    db = MitoVariantDatabase()
    if args.rcs:
        get_or_create_mmi(args.rcs)
    if args.add:
        if not all([args.bam, args.acq, args.pid, args.rcs, args.outdir, args.db]):
            print("Error: -add requires --bam, --acq, --pid, --outdir, --db, and --rcs arguments.")
            sys.exit(1)
            
        try:
            db.load(args.db)
        except FileNotFoundError:
            print(f"no db loaded, doing default")

        print(f"adding {args.bam}")
        if is_unaligned_bam(args.bam):
            vcf_file, cov = make_vcf_unaligned_bam(args.bam, args.rcs, args.outdir)
        else:
            vcf_file = make_vcf(args.bam, args.rcs, args.outdir)
            cov = get_coverage(args.bam)
        
        if vcf_file:
            db.add_sample(vcf_file, args.acq, args.pid, cov)
            db.save(args.db)
            print(f"Added sample {args.acq} to database")
        remove_files(args.outdir)

            
    elif args.remove:
        if not all([args.acq, args.db]):
            print("Error: -remove requires the --acq and --db")
            sys.exit(1)
            
        try:
            db.load(args.db)
        except FileNotFoundError:
            print(f"Error: no database to remove from, tried {args.db}")
            sys.exit(1)
            
        print(f"Removing sample {args.acq}...")
        success = db.remove_sample(args.acq)
        if success:
            db.save(args.db)
            print(f"updated sample saved to {args.db}.")
        else:
            print("sample not removed")

    elif args.create:
        if not all([args.csv, args.rcs, args.outdir, args.db]):
            print("Error: -create requires --csv --db --outdir and --rcs arguments.")
            sys.exit(1)
            
        print("Creating database from CSV, give it a few mins")
        os.makedirs(args.outdir, exist_ok=True)
        db.create_database_from_csv(args.rcs, args.outdir, args.csv, args.db)
        print(f"Database saved to {args.db}")
        remove_files(args.outdir)


    elif args.compare:
        if not all([args.bam or args.filter_sample, args.rcs, args.outdir, args.db, args.threshold]):
            print("Error: -compare requires --bam or --filter_sample, --rcs, --outdir, --threshold, --db arguments.")
            sys.exit(1)
            
        try:
            db.load(args.db)
        except FileNotFoundError:
            print(f"Error: no database {args.db}")
            sys.exit(1)
        
        if args.filter_sample:
            #THIS WILL JUST GIVE FOR ALL THE MATCHES. It's honestly impossible to check for only the ones that this sample matches without doing checkall, so why
            #not just do checkall after every add, and if you want to look for a specific sample just search for the name
            #This will give you enough, a check on the matches and the supoosed ID given if it's already in the database (which it needs to be)
            print(f"comparing {args.filter_sample}")
            db.make_igv_output_putative_compare(args.threshold, args.outdir, args.filter_sample)
        
        elif args.bam:
            print(f"comparing {args.bam}")
            if is_unaligned_bam(args.bam):
                vcf_file, cov = make_vcf_unaligned_bam(args.bam, args.rcs, args.outdir)
            else:
                vcf_file = make_vcf(args.bam, args.rcs, args.outdir)
                cov = get_coverage(args.bam)
                
            if vcf_file:
                print("getting matches")
                matches = db.compare_sample(vcf_file, threshold=args.threshold)
                
                # Use placeholder IDs for the external sample evaluation
                db.make_output(matches, db, args.bam, "", args.outdir)
                
                print("visualization")
                
                #print(args.compare)
                outfile = args.outdir + "/" + args.visualize
                # print(args.compare_sample)
                
                #make putative does this now
                if args.bam:
                    db.visualize_sample(vcf_file, outfile = outfile, coverage=cov, threshold=args.threshold, compare = args.compare_sample)
                
                
            
            print(f"visualization output to {outfile}")
        remove_files(args.outdir)


    elif args.checkall:
        if not all([args.outdir, args.db]):
            print("Error: -checkall requires --outdir and --db arguments.")
            sys.exit(1)
        try:
            db.load(args.db)
        except FileNotFoundError:
            print(f"Error: Could not find existing database at {args.db}")
            sys.exit(1)
            
        print("comparing all samples in database")
        db.check_all_database(args.outdir, args.hist, args.log_file, args.threshold, args.vis_putative)
        print("histogram generated and matches logged")
        remove_files(args.outdir)


if __name__ == "__main__":
    main()
