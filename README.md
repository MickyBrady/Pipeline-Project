Computational Biology Pipeline Project

#Pipeline for Transcriptome Analysis

#Step 1: Retrieving Transcriptomes from SRA and Converting to FASTQ Files

We first retrieve the transcriptomes from SRA (Sequence Read Archive) for the patient donors to begin the analysis pipeline. This step involves downloading the SRA files and converting them to paired-end FASTQ files using `fastq-dump`.

# Donors and Time Points
- **Donor 1 (2dpi)**: SRX2896360
- **Donor 1 (6dpi)**: SRX2896363
- **Donor 3 (2dpi)**: SRX2896374
- **Donor 3 (6dpi)**: SRX2896375

# Commands Used:
We used the `prefetch` and `fastq-dump` tools to retrieve and convert the SRA files into FASTQ format.

# Step 1: Fetching the SRA Files
prefetch SRX2896360
prefetch SRX2896363
prefetch SRX2896374
prefetch SRX2896375
# Converting SRA files to fwd and reverse read fastq files using fastq-dump 
fastq-dump --split-files SRX2896360
fastq-dump --split-files SRX2896363
fastq-dump --split-files SRX2896374
fastq-dump --split-files SRX2896375

Step 1: Retrieve Transcriptomes from SRA
Overview
For this project, you will retrieve the transcriptomes of two patient donors at 2- and 6-days post-infection (dpi) from the SRA (Sequence Read Archive). These datasets are publicly available on NCBI and will be converted into paired-end FASTQ files for downstream analysis. The four relevant datasets are:

Donor 1 (2dpi): SRX2896360
Donor 1 (6dpi): SRX2896363
Donor 3 (2dpi): SRX2896374
Donor 3 (6dpi): SRX2896375
These files can be retrieved using either the SRA toolkit or wget.


