import sys
import argparse
from collections import defaultdict
import pickle
from Bio import SeqIO

def find_cpg_sites(seq: str) -> list[int]:
    """Return 0-based positions of all CpG sites in a sequence."""
    cpg_sites = []
    for i in range(len(seq) - 1):
        if seq[i:i+2].upper() == "CG":
            cpg_sites.append(i)
    return cpg_sites

def extract_cpg_from_fasta(fasta_path: str) -> dict:
    """Extract CpG site positions from a FASTA file."""
    cpg_dict = defaultdict(list)
    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = record.id
        if "_" in chrom:
            continue
        seq = str(record.seq)
        positions = find_cpg_sites(seq)
        cpg_dict[chrom] = positions
        print(f"Found {len(positions)} CpG sites in {chrom}")
    return dict(cpg_dict)

def save_cpg_sites(cpg_sites: dict, output_path: str):
    """Save the CpG site dictionary to a pickle file."""
    with open(output_path, "wb") as f:
        pickle.dump(cpg_sites, f)
    print(f"Saved CpG site positions to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Extract CpG sites from a FASTA reference.")
    parser.add_argument("fasta", help="Reference genome in FASTA format")
    parser.add_argument("-o", "--output", default="cpg_sites.pkl", help="Output file (default: cpg_sites.pkl)")
    args = parser.parse_args()

    cpg_sites = extract_cpg_from_fasta(args.fasta)
    save_cpg_sites(cpg_sites, args.output)

if __name__ == "__main__":
    main()
