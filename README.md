# Cell-Ranger
Cell Ranger - From Algorithm to working tutorial

Cell Ranger is used to process data generated from 10x Genomics single-cell experiments.

  •	The mkfastq command converts raw BCL files (from Illumina sequencers) into FASTQ files.
  
  •	The count command processes these FASTQ files to generate a gene–barcode expression matrix.
  
  •	The aggr (aggregate) command combines results from multiple count runs, enabling joint analysis and normalization across samples.
________________________________________


# Cell Ranger algorithm:

FastQ (R1,R2,I1/I2)

      │
      ▼
      
Barcode & UMI extraction

      │
      ▼
      
Read alignment to reference (STAR)

      │
      ▼
      
UMI deduplication & gene assignment

      │
      ▼
      
Cell calling (distinguish real cells from empty droplets)

      │
      ▼
      
Gene x Cell matrix + summary metrics


1. Barcode & UMI extraction (from Read 1)
   
•	The R1 read in 10x libraries contains:

o	Cell barcode (usually 16 bp)

o	UMI (Unique Molecular Identifier) (usually 10–12 bp)

•	These sequences are extracted before any alignment happens.

•	Cell Ranger:

o	Reads R1 to get the barcode and UMI.

o	Reads R2 to get the cDNA (the actual transcript sequence).

o	Links each R2 read to its barcode and UMI.

 Why first?
Because R1 (barcode+UMI) is not biological transcript — it’s an artificial sequence added during library prep. It doesn’t align to the genome.
You must extract it first so you know which cell and molecule each R2 belongs to before mapping.

2.  Alignment (of R2 reads)
   
•	Only the R2 sequences are aligned to the reference genome/transcriptome.

•	Cell Ranger uses STAR aligner internally.

•	After alignment, it knows:

o	Which gene each R2 read maps to.

o	Which barcode (cell) and UMI (molecule) that read belongs to.

•	Reads that map ambiguously can be filtered or assigned probabilistically.

3.  Post-alignment processing
   
After alignment, Cell Ranger combines the information:

•	Barcode → identifies the cell.

•	UMI → identifies the molecule (unique transcript copy).

•	Gene → identified from alignment.

Then it does:

•	Barcode correction:

o	10x barcodes can have sequencing errors.

o	Cell Ranger compares barcodes to the whitelist (known 10x barcodes) and corrects single-base errors.

•	UMI deduplication:

o	Reads with the same UMI and same gene are counted only once.

o	Reduces PCR duplication bias.

Result: A table of unique molecules per gene per cell.

4.  Cell Calling / Filtering
   
•	Purpose: Distinguish real cells from empty droplets.

•	Algorithm:

o	Calculate total UMI counts per barcode.

o	Fit a distribution of counts (e.g., a knee plot).

o	Barcodes with low counts are considered empty droplets; high counts are real cells.

Output: List of high-confidence cells and their counts.

________________________________________

# CELL RANGER WORKFLOW

1.  Installation
   
a. Download Cell Ranger

From the 10x Genomics website:
 https://www.10xgenomics.com/support/software/cell-ranger
 
bash:
```
cellranger-9.0.1.tar.gz
```
b. Extract the tar file
```
tar -xzvf cellranger-9.0.1.tar.gz
```
This creates a folder : cellranger-9.0.1/

c. Call it using the full path, e.g.:
```
~/cellranger-9.0.1/bin/cellranger
```
2️. Prepare Input Data

a. Download the reference genome

Get the prebuilt 10x reference for your species, e.g.:
•	Human (GRCh38):
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

Bash:
```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
```
Extract it:
```
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```
You’ll get a folder like:

refdata-gex-GRCh38-2020-A/

3. Download  dataset
   
You can use datasets from the 10x website.

bash:
```
wget https://cf.10xgenomics.com/samples/cell-exp/6.0.0/1k_PBMCs_TotalSeq_B_3p/1k_PBMCs_TotalSeq_B_3p_fastqs.tar
```
Extract it:
```
tar -xzvf 1k_PBMCs_3p_v3_fastqs.tar.gz
```
You’ll get a folder with .fastq.gz files (R1, R2, I1, I2).

4. Run cellranger count
   
Command syntax:
```
home/bioinformatics/cellranger-9.0.1/bin/cellranger count --id run1output --transcriptome ./refdata-gex-GRCh38-2024-A/ --fastqs ./1k_PBMCs_TotalSeq_B_3p_fastqs/1k_PBMCs_TotalSeq_B_3p_gex_fastqs/ --sample 1k_PBMCs_TotalSeq_B_3p_gex --create-bam=false
```
Parameters :

•	--id → name for output folder

•	--transcriptome → path to reference genome folder

•	--fastqs → path to folder containing FASTQ files

•	--sample → prefix of sample (matches the beginning of your FASTQ filenames)

5. Output
   
Run1output file contains output which include

Gene-barcode matrix (primary output for downstream analysis, e.g., Seurat)

BAM file (optional, if --create-bam=true)

Summary metrics:

o	Number of cells

o	Median genes per cell

o	Sequencing saturation

o	Fraction of reads mapped


