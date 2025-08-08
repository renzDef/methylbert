# Information for using the MethylBERT model

### Scripts

- `convert_cpg_sites.py`: Creates the cpg_index_to_pos from the .pkl file, see comment at additional files.
- `dmr_calling.py`: Convertes the information from the bigwig data to DMR data and saves to .csv file. MethylBERT people used some R tool, but that did not work for me, maybe I just did not understand R. Usage:
<pre><code>
python dmr_calling.py ../data/bigwig ../data/groups.txt ../data/reference/dmr.csv
</code></pre>
- `extract_cell_types.py`: Checks the directory with all pad files and extracts the names of the cells.
- `extract_methylation_sites.py`: Extracts all CpG sites from the reference genome. Usage:
<pre><code>
python extract_methylation_sites.py ../data/reference/hg38.fa
</code></pre>
- `filter_dmrs.py`: Used to change the `dmr.csv` in a way so that the Methylseq Simulation works with that. (Not used right now, but wanted to generate the bulk data with that, maybe will change it later, but did not work as expected.)
- `generate_bulk_sample.py`: Generates bulk sample data by combinig reads from the different BAM files specific for the cell types. Currently uses only the first to classes.
- `pat_to_sam.py`: Creates reads from the information of the pad file. Usage:
<pre><code>
python pat_to_sam.py ../data/pat/name_of_pat_file.pat ../data/reference/cpg_sites.pkl ../data/reference/hg38.fa ../data/bam_for_fine_tuning/name_of_bam_file.bam
</code></pre>
- `process_pat_files.sh`: Script for autmatically converting all pat files in the pat directory to bam files. (Skips extisting files). Usage:
<pre><code> ./process_pat_files.sh ../data/pat ../data/bam_for_fine_tuning</code></pre>



### Data

- bam_for_classification: Created with `generate_bulk_sample.py`
- bam_for_fine_tuning: Created with `proces_pat_files.sh`
- bigwig: downloaded, used for extracting DMRs in `dmr_calling.py`
- pat: downloaded, used to generate the bam files for fine-tuning.
- reference: `hg38.fa` (downloaded), `cpg_sites.pkl` (extracted from reference genome), `dmr.csv` (created with `dmr_calling.py`), `dmr_filtered.csv` (version with only longer reads and less rows)
- groups.txt: list of the classes the biwig data belongs to, used for DMR calling.

### Additional Files

- cell_types.tsv: List of all names of cell types and the class
- cpg_index_to_pos.tsv: I forgot what I used that for, probably in the process of getting the DMR data. Will check and remove if useless now.

## Conda env

to recreate the conda environment, run the following code in the `methylbert` folder:

```
conda env create -f env_config.yml -n methylbert

```

After this you should be able to run all python files.
