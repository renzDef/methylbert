''' DMR data is calculated from the BigWig data files.
The following should be provided:
- path to bigwig data
- file with the classes of the data in the order of processing (groups.txt)
- output csv file

'''

#!/usr/bin/env python3

import argparse
import pyBigWig
import pandas as pd
from collections import defaultdict
import subprocess
import tempfile
import os

def read_groups_file(groups_path):
    with open(groups_path) as f:
        return [int(x) for x in f.read().strip().split()]

def extract_common_cpgs(bigwig_files, max_sites=1_000_000):
    """Extract common CpG positions across all BigWig files"""
    site_sets = []
    for bw_path in bigwig_files:
        bw = pyBigWig.open(bw_path)
        sites = set()
        print(f"  Reading {os.path.basename(bw_path)}...")
        for chrom in bw.chroms().keys():
            values = bw.intervals(chrom)
            if values:
                for start, end, val in values:
                    sites.add((chrom, start))
                    if len(sites) >= max_sites:
                        break
        bw.close()
        print(f"    -> {len(sites)} sites extracted.")
        site_sets.append(sites)

    if not site_sets:
        raise ValueError("No CpG sites found in any BigWig file.")

    common_sites = set.intersection(*site_sets)
    if not common_sites:
        raise ValueError("No common CpG sites found across all BigWig files.")

    print(f"  => {len(common_sites)} common CpG sites found across all samples.")
    return sorted(common_sites)

def extract_methylation_matrix(bigwig_files, sites, groups):
    """Create a methylation matrix from BigWigs at the given sites"""
    methylation_dict = defaultdict(list)
    sample_names = [os.path.splitext(os.path.basename(f))[0] for f in bigwig_files]
    sample_names_with_group = ["g" + str(groups[i]) + "_" + sample_names[i] for i in range(len(groups))]

    for bw_path in bigwig_files:
        bw = pyBigWig.open(bw_path)
        for chrom, start in sites:
            val = bw.values(chrom, start, start+1)[0]
            val = 0.0 if val != val else val
            methylation_dict[(chrom, start)].append(val)
        bw.close()

    rows = []
    for (chrom, start), values in methylation_dict.items():
        rows.append([chrom, start] + values)
    columns = ["chr", "pos"] + sample_names_with_group
    df = pd.DataFrame(rows, columns=columns)
    return df

def run_metilene(matrix_path, groups_path, output_path):
    try:
        with open(output_path, "w") as out:
            subprocess.run(["metilene", "-a" , "g0", "-b", "g1", matrix_path], stdout=out, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Metilene failed: {e}")

def format_dmrs_for_methylbert(dmrs_path, output_csv_path):
    df = pd.read_csv(dmrs_path, sep='\t', header=None)
    df.columns = ["chr", "start", "end", "nCG", "meanMethy1", "meanMethy2", "diff.Methy", "p", "q", "areaStat"]
    df["length"] = df["end"] - df["start"]
    df["ctype"] = "group_comparison"
    df = df[["chr", "start", "end", "length", "nCG", "meanMethy1", "meanMethy2", "diff.Methy", "areaStat", "ctype"]]
    df.to_csv(output_csv_path, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Call DMRs from BigWig files using metilene for MethylBERT.")
    parser.add_argument("bigwig_dir", help="Directory with BigWig files")
    parser.add_argument("groups_file", help="groups.txt file for metilene")
    parser.add_argument("output_file", help="Output DMR file for MethylBERT")
    args = parser.parse_args()

    bigwig_files = [args.bigwig_dir + f for f in os.listdir(args.bigwig_dir) if os.path.isfile(os.path.join(args.bigwig_dir, f))]
    groups = read_groups_file(args.groups_file)

    if len(groups) != len(bigwig_files):
        raise ValueError(f"groups.txt must have {len(bigwig_files)} values, but has {len(groups)}.")

    print("[1/4] Extracting common CpG sites...")
    common_sites = extract_common_cpgs(bigwig_files)

    print("[2/4] Building methylation matrix...")
    matrix_df = extract_methylation_matrix(bigwig_files, common_sites, groups)
    print(matrix_df.columns)

    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".tsv") as tmp_matrix:
        matrix_df.to_csv(tmp_matrix.name, sep='\t', index=False, header=True)
        tmp_matrix_path = tmp_matrix.name

    print("[3/4] Running metilene...")
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".tsv") as tmp_dmrs:
        run_metilene(tmp_matrix_path, args.groups_file, tmp_dmrs.name)

        print("[4/4] Formatting DMRs for MethylBERT...")
        format_dmrs_for_methylbert(tmp_dmrs.name, args.output_file)

    print(f"✅ DMR file written to: {args.output_file}")

if __name__ == "__main__":
    main()
