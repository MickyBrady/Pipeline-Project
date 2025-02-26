# Computational Biology Pipeline Project

This project aims to develop a Python wrapper script to automate the execution of various bioinformatics tools for analyzing HCMV transcriptomes.

# List of required Tools 
Python3 

Biopython

Kallisto

Bowtie2

SPAdes

BLAST + 

R and Sleuth package

# Wrapper Python Script

to run this script, please clone this GitHub repository with: 

```bash
git clone https://github.com/MickyBrady/Pipeline-Project.git
```

This repository contains a sample_data set of transcriptomes from two patient donors. 

Using wget this sample took the first 10,000 reads from  SRA and converted to paired-end fastq files.

```
wget
Donor 1 (2dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
Donor 1 (6dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
Donor 3 (2dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
Donor 3 (6dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
```






