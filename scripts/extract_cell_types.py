''' script to extract names of the cells which cell types are in the data. '''

import os
import re
import csv

def extract_unique_middle_parts(directory, output_file):
    unique_names = set()

    # Pattern to match the filename format and extract the middle part
    pattern = re.compile(r'^GSM\d+_([^-]+)-')

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            unique_names.add(match.group(1))

    # Write the unique names to a tab-separated CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for name in sorted(unique_names):
            writer.writerow([name])

    print(f"Extracted {len(unique_names)} unique names to '{output_file}'.")

# Example usage
# Provide your directory and output file path here
extract_unique_middle_parts('data/pat/', 'cell_types.tsv')
