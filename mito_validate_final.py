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


def make_vcf(bam_file, rcs_file, output_dir, mt_chrom = "chrM", level = ".9", post_filter = True):
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
    mt_bam = os.path.join(output_dir, "temp.bam")
    vcf_output = os.path.join(output_dir, "temp_result.vcf.gz")
    try:
        # remove old VCF if it exists so Mutserve doesn't complain about overwriting
        if os.path.exists(vcf_output):
            os.remove(vcf_output)

        # extract Mitochondrial Reads to the temp BAM
        subprocess.run([
            "samtools", "view", "-b",
            "-F", "2048",       # exclude Supplementary
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
            "--level", level, #Should adjust this number
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
                    # Iterate over samples (usually just one in Mutserve output)
                    # logic: if ANY sample has AF >= 0.9, keep the record
                    keep_record = False
                    
                    for sample in record.samples.values():
                        # Get AF, default to 0.0 if missing
                        af = sample.get('AF', 0.0)
                        
                        # Pysam often returns AF as a tuple/list
                        if isinstance(af, (list, tuple)):
                            af = af[0]
                        
                        if float(af) >= 0.9:
                            keep_record = True
                            break
                    
                    if keep_record:
                        vcf_out.write(record)
                
                vcf_in.close()
                vcf_out.close()
                
                #overwrite the original output with the filtered one
                os.replace(filtered_vcf, vcf_output)            
            
        
        
        #remove bam
        if os.path.exists(mt_bam):
            os.remove(mt_bam)
        if os.path.exists(mt_bam + ".bai"):
            os.remove(mt_bam + ".bai")
        

        return vcf_output

    except subprocess.CalledProcessError as e:
        print(f"Failed to process {bam_file}: {e}")
        return None

def get_coverage(mt_bam, mt_chrom = "chrM"):
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
        # e.stderr contains the error message from samtools
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

    def add_sample(self, vcf_path, sample_id, bank_id, coverage=None):
        """
        Adds a new sample to the existing sparse matrix database.
        
        Inputs:
            vcf_path (str): Path to the VCF file to be added.
            sample_id (str): Unique identifier for the acquisition/sample.
            bank_id (str): Identifier for the biological source or bank.
            coverage (float): Mean sequencing depth of the sample.
            
        Outputs:
            None. Modifies self.db_matrix, self.sample_ids, and self.cash_sum in place.
        """
        if ([sample_id, str(bank_id), coverage] in self.sample_ids):
            print("sample already in database, not actually adding")
        else:
            new_vector = self.vcf_to_vector(vcf_path)
            # add the new column to existing DB
            self.db_matrix = hstack([self.db_matrix, new_vector])
            self.sample_ids.append((sample_id, str(bank_id), coverage))
            
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
        
        print(f"Successfully removed sample: {removed_metadata[0]} (Bank: {removed_metadata[1]})")
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

        sample_vec = self.vcf_to_vector(vcf_path)
        
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
    
    
    def visualize_sample(self, vcf_path, outfile, coverage = 100.0, threshold=0.8, compare = None): #set to 100 so it works so the default isn't low coverage
        """
        #get the database where the columns are the ones that match
        
        #make a plot like igv, where the the length is 16k bases, and each snp in sample_vec is shown.
        #I want each of the matching samples to be shown above with the matching snps highlighted green and non matching ones red
        """
        
        """
        Generates an IGV-style plot comparing a target sample's variants against matching database samples.
        
        Inputs:
            vcf_path (str): Path to the target VCF file.
            coverage (float): Coverage of the target sample; influences plotting color (default: 100.0).
            threshold (float): Minimum similarity threshold to retrieve matches for plotting (default: 0.8).
            compare (str): Specific bank_id to force a comparison with, overriding the threshold (default: None).
            
        Outputs:
            None. Renders a matplotlib plot to the screen.
        """
        

        # 1. Parse the input sample
        sample_vec = self.vcf_to_vector(vcf_path)
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
        y_labels = ["Target Input"]

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
        ax.set_title(f"Variant Comparison (Threshold >= {threshold})")
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
    
    def make_output(self, matches, mito_db, acq_id, bank_id, output_dir, file_name = "matches_log.txt"):
        """
        Appends matching samples to a text file.
        
        Inputs:
            matches (numpy.ndarray): Boolean array indicating matched samples.
            mito_db (MitoVariantDatabase): The database instance containing sample metadata.
            acq_id (str): Acquisition ID of the evaluated sample.
            bank_id (str): Bank ID of the evaluated sample.
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
            log_message = f"Sample {acq_id}, {bank_id} matches existing samples: {matching_names}\n"
                        
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
            sample_bam = acq_folder + f"/{acq_id}.sorted.bam"
            bank_id = str(row.name)
            
            vcf_file = make_vcf(sample_bam, rcs_file, output_dir)
            
            coverage = get_coverage(sample_bam)
            #check matches
            #matches = mito_db.compare_sample(vcf_file, threshold=THRESHOLD) 
            
            #make_output(matches, mito_db, acq_id, bank_id, output_dir)

            # add to Database
            self.add_sample(vcf_file, sample_id=acq_id, bank_id = bank_id, coverage = coverage)

        #save Final Database
        self.save(save_file)


    def check_all_database(self, output_dir, histogram_save, log_file, threshold):
        """
        Calculates pairwise similarity for all samples currently in the database 
        and plots a histogram of the scores. Excludes self-comparisons.
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
            for i, j in zip(upper_tri_idx[0], upper_tri_idx[1]):
                coverage_warning = ""
                cov_end = ""
                if self.sample_ids[i][2] < 10:
                    coverage_warning += f"Warning: coverage on sample {self.sample_ids[i][0]} = {self.sample_ids[i][2]} < 10\n"
                    cov_end = "\n"
                if  self.sample_ids[j][2] < 10:
                    coverage_warning += f"Warning: coverage on sample {self.sample_ids[j][0]} = {self.sample_ids[j][2]} < 10\n"
                    cov_end = "\n"
                if scores_matrix[i, j] >= threshold:
                    if self.sample_ids[i][1] != self.sample_ids[j][1]:
                        # sample_ids elements are tuples like (sample_id, bank_id, coverage)
                        s1 = f"{self.sample_ids[i][0]} ({self.sample_ids[i][1]})"
                        s2 = f"{self.sample_ids[j][0]} ({self.sample_ids[j][1]})"
                        f.write(f"{cov_end}{coverage_warning}Match above threshold for different Bank IDs ({scores_matrix[i, j]:.2f}): {s1} <--> {s2}{cov_end}\n")
                    elif scores_matrix[i, j] < 1:
                        s1 = f"{self.sample_ids[i][0]} ({self.sample_ids[i][1]})"
                        s2 = f"{self.sample_ids[j][0]} ({self.sample_ids[j][1]})"
                        f.write(f"{cov_end}{coverage_warning}Match above threshold but below 1 same Bank ID ({scores_matrix[i, j]:.2f}): {s1} <--> {s2}{cov_end}\n")
                elif self.sample_ids[i][1] == self.sample_ids[j][1]:
                    s1 = f"{self.sample_ids[i][0]} ({self.sample_ids[i][1]})"
                    s2 = f"{self.sample_ids[j][0]} ({self.sample_ids[j][1]})"
                    f.write(f"{cov_end}{coverage_warning}Match below threshold from same Bank ID ({scores_matrix[i, j]:.2f}): {s1} <--> {s2}{cov_end}\n")

        # plot
        passing_scores = all_scores[all_scores >= threshold]
        failing_scores = all_scores[all_scores < threshold]
            
        plt.figure(figsize=(10, 6))
            
        plt.hist([passing_scores, failing_scores], 
                bins=50, 
                stacked=True, 
                color=['green', 'red'], 
                label=['Above Threshold', 'Below Threshold'],
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
    parser.add_argument('--rcs', type=str, help='Path to reference FASTA (Required for -add, -create, -compare)')
    parser.add_argument('--threshold', type=float, default=0.8, help='Similarity threshold (default: 0.8)')
    
    # main helpers
    parser.add_argument('--csv', type=str, help='Path to the CSV file (Required for -create)')
    parser.add_argument('--bam', type=str, help='Path to the BAM file (Required for -add and -compare)')
    parser.add_argument('--acq', type=str, help='Acquisition ID (Required for -add and -remove) This is the individual sample ID')
    parser.add_argument('--bank', type=str, help='Bank ID (Required for -add) This is the group ID, which can contain multiple acquisition IDs')
    parser.add_argument('--hist', type=str, default='score_histogram.png', help='Histogram output filename (Used with -checkall)')
    parser.add_argument('--log_file', type=str, default='pairwise_matches_log.txt', help='Name of log file (Used with -checkall)')
    parser.add_argument('--visualize', type=str, default='visualization.png', help='Visualize file matches after comparing (Optional for -compare)')
    parser.add_argument('--compare_sample', type=str, default = "", help='A specific bank ID to compare a sample to (not just the ones it matches) (Optional for -compare)')

    args = parser.parse_args()
    db = MitoVariantDatabase()

    if args.add:
        if not all([args.bam, args.acq, args.bank, args.rcs, args.outdir, args.db]):
            print("Error: -add requires --bam, --acq, --bank, --outdir, --db, and --rcs arguments.")
            sys.exit(1)
            
        try:
            db.load(args.db)
        except FileNotFoundError:
            print(f"no db loaded, doing default")

        print(f"adding {args.bam}")
        vcf_file = make_vcf(args.bam, args.rcs, args.outdir)
        if vcf_file:
            cov = get_coverage(args.bam)
            db.add_sample(vcf_file, args.acq, args.bank, cov)
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
        if not all([args.bam, args.rcs, args.outdir, args.db]):
            print("Error: -compare requires --bam, --rcs, --outdir, --db arguments.")
            sys.exit(1)
            
        try:
            db.load(args.db)
        except FileNotFoundError:
            print(f"Error: no database {args.db}")
            sys.exit(1)
        
        print(f"comparing {args.bam}")
        vcf_file = make_vcf(args.bam, args.rcs, args.outdir)
        if vcf_file:
            print("getting matches")
            matches = db.compare_sample(vcf_file, threshold=args.threshold)
            
            # Use placeholder IDs for the external sample evaluation
            db.make_output(matches, db, args.bam, "", args.outdir)
            
            print("visualization")
            cov = get_coverage(args.bam)
            #print(args.compare)
            outfile = args.outdir + "/" + args.visualize
            # print(args.compare_sample)
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
        db.check_all_database(args.outdir, args.hist, args.log_file, args.threshold)
        print("histogram generated and matches logged")
        remove_files(args.outdir)


if __name__ == "__main__":
    main()
