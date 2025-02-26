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

# Clone Repository 

to run this script, please clone this GitHub repository with: 

```bash
git clone https://github.com/MickyBrady/Pipeline-Project.git
```

This repository includes a smaller sample_data subset of transcriptomes from two patient donors.
You can download the full dataset from SRA using the following wget commands. 
The first 10,000 reads from each sample were converted to paired-end fastq files for quick testing.
```
wget
Donor 1 (2dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
Donor 1 (6dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
Donor 3 (2dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
Donor 3 (6dpi): https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
```
# Running the Wrapper.py 

Once you've cloned the repository and downloaded the data, 
you can run the pipeline using the wrapper.py script.

To execute the pipeline, navigate to the directory and use the following command:

 ```bashh
cd Pipeline-Project
```

```bashh
python3 wrapper.py
```






