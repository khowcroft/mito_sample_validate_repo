Activate environment mito_validate_vcfs_matrix

Note: for the database prefix --db: if the file name is mitochondria_sample_validation/final/output/add_to_none_test_meta.json, the prefix is add_to_none_test

## Create Database:
Takes a CSV and creates a database of mitochondrial SNPs
CSV format (important is Bank ID, acquisitions, acq_id)
Bank ID,Group,Brain Region,acquisitions,acq_id,pass-reads,n50,LINE_single_end,LINE_double_end,LINE_none,AluYa5_single_end,AluYa5_double_end,AluYa5_none,AluYb8_single_end,AluYb8_double_end,AluYb8_none,LINE_enrichment,AluYa5_enrichment AluYb8_enrichment,Known_polymorphic_LINE-1_captured,Known_polymorphic_AluY_captured
196,Control,Frontal Cortex,/path_to/acquisitions,acq_id,"236,900","7,571",35149,12005,177792,66625,11899,146422,64887,10661,149398,20.96%,34.91%,33.58%,187,679

-create requires --csv --db --outdir and --rcs arguments

ex.
    python /path_to/mito_validate_final.py -create --rcs /path_to/rCRS.fasta --csv /path_to/experiments.csv --db /path_to/output/test_database_prefix --outdir /path_to/output

## Add sample:
Adds a sample to a new database or an existing database
-add requires --bam, --acq, --bank, --outdir, --db, and --rcs arguments

ex.
    python /path_to/mito_validate_final.py -add --bam /path_to/sample.sorted.bam --acq sample1_id --bank person_id --outdir /path_to/output --db /path_to/output/db_to_add_prefix --rcs /path_to/rCRS.fasta

## Remove Sample:
removes a sample from the database
-remove requires the --acq and --db

ex.
    python /path_to/mito_validate_final.py -remove --db /path_to/test_database_prefix --acq sample_id_to_remove

## Check all samples in a database:
compares all the samples in the database and gives a histogram of all the comparison scores and an output file describing which samples match eachother
Specifically, it will give:
1. Samples where the bank IDs match but a score of 1 (identical) isn't achieved but is still above the threshold]
2. Samples where the bank IDs match and the score is below the threshold
3. Samples where the bank IDs dont match but the score is above the threshold
-checkall requires --outdir and --db arguments
takes an optional --hist  and or --threshold and or --log_file 
    
ex.
    python /path_to/mito_validate_final.py -checkall --outdir /path_to/output --db /path_to/output/database_prefix
or with optional args:
    python /path_to/mito_validate_final.py -checkall --outdir /path_to/output --hist hist2.png --threshold .85 --log_file log_file.txt --db /path_to/output/database_prefix


## Compare a sample to the database without adding:
Compares the given sample to all of the samples in the database and outputs a line to matches.txt in outdir (accumulates) as well as the vizualization compared to the matches or the bank ID of interest
-compare requires --bam, --rcs, --outdir, --db arguments
optionally can use --visualize to name the visualization output and --compare_sample to name a specific bank_id to compare to

ex.
    python /path_to/mito_validate_final.py -compare --bam /path_to/sample.sorted.bam --rcs /path_to/rCRS.fasta --outdir /path_to/output --db /path_to/test_database_prefix --visualize vis_file.png --compare_sample bankID_name




## args
```
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
```
