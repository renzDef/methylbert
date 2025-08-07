''' This generates the data needed for testing the fine-tuned model by combining reads from generated bam files.'''

import pysam
import random

def select_random_reads(bam_file, n_reads):
    """
    Randomly select n_reads from a BAM file.

    Parameters:
    bam_file (str): The path to the BAM file.
    n_reads (int): The number of reads to randomly select.

    Returns:
    list: List of selected reads.
    """
    # Open the BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")

    # Get all reads
    all_reads = [read for read in samfile.fetch()]

    # Randomly select n_reads
    selected_reads = random.sample(all_reads, min(n_reads, len(all_reads)))

    return selected_reads

def combine_bam_files(bam_file_0, bam_file_1, n_reads_0, n_reads_1, output_file):
    """
    Combine random reads from two BAM files into one.

    Parameters:
    bam_file_0 (str): The path to the BAM file for cell type 0.
    bam_file_1 (str): The path to the BAM file for cell type 1.
    n_reads_0 (int): The number of reads to select from bam_file_0.
    n_reads_1 (int): The number of reads to select from bam_file_1.
    output_file (str): The path to save the combined BAM file.
    """
    # Select random reads from both BAM files
    selected_reads_0 = select_random_reads(bam_file_0, n_reads_0)
    selected_reads_1 = select_random_reads(bam_file_1, n_reads_1)

    # Open the output BAM file for writing
    with pysam.AlignmentFile(output_file, "wb", header=pysam.AlignmentFile(bam_file_0, "rb").header) as out_bam:
        # Write the selected reads from both cell types into the output BAM file
        for read in selected_reads_0 + selected_reads_1:
            out_bam.write(read)

# Example usage
bam_file_0 = "../data/bam_for_fine_tuning/GSM5652176_Adipocytes-Z000000T7.bam"
bam_file_1 = "../data/bam_for_fine_tuning/GSM5652179_Aorta-Endothel-Z00000422.bam"
n_reads_0 = 1000  # Select 1000 reads from cell type 0
n_reads_1 = 1000  # Select 1000 reads from cell type 1
output_file = "../data/bam_for_classification/combined_bulk_sample.bam"  # Combined BAM file with reads from both cell types

combine_bam_files(bam_file_0, bam_file_1, n_reads_0, n_reads_1, output_file)
