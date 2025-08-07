# ''' script to convert cpg cites from .pkl format to .csv format. '''


import pickle
import pandas as pd

# Load your existing CpG sites file
with open("../data/reference/cpg_sites.pkl", "rb") as f:
    cpg_dict = pickle.load(f)  # Structure: { "chr1": [0, 1, 5, ...], ... }

# Convert to flat table with global 1-based index and 1-based positions
records = []
index = 1
for chrom in sorted(cpg_dict.keys()):
    for pos in cpg_dict[chrom]:
        records.append((index, chrom, pos + 1))  # convert to 1-based genome position
        index += 1

# Create DataFrame
df = pd.DataFrame(records, columns=["index", "chr", "pos"])

# Save as TSV for R
df.to_csv("cpg_index_to_pos.tsv", sep="\t", index=False)
print(f"Saved {len(df)} CpG sites to cpg_index_to_pos.tsv")
