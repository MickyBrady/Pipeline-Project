import os
import subprocess

# Define paths
pipeline_dir = "/home/2025/mbrady9/PipelineProject_Michaela_Brady"
sample_data_dir = os.path.join(pipeline_dir, "sample_data")
log_file_path = os.path.join(pipeline_dir, "PipelineProject.log")
fasta_filename = os.path.join(pipeline_dir, "HCMV_cds.fasta")
index_prefix = os.path.join(pipeline_dir, "HCMV_index")

# Function to write the number of reads before and after Bowtie2 mapping to the log file
def write_bowtie2_log(log_file, donor, condition, reads_before, reads_after):
    with open(log_file, 'a') as log:
        log.write(f"{donor} ({condition}) had {reads_before} read pairs before Bowtie2 filtering and {reads_after} read pairs after.\n")

# Create Bowtie2 genome index
bowtie2_index_command = f"bowtie2-build {fasta_filename} {index_prefix}"
subprocess.run(bowtie2_index_command, shell=True, check=True)

# Define sample FASTQ files and conditions
samples = {
    "SRR5660030": (os.path.join(sample_data_dir, "sample_SRR5660030_1.fastq"), os.path.join(sample_data_dir, "sample_SRR5660030_2.fastq"), "Donor 1", "2dpi"),
    "SRR5660033": (os.path.join(sample_data_dir, "sample_SRR5660033_1.fastq"), os.path.join(sample_data_dir, "sample_SRR5660033_2.fastq"), "Donor 1", "6dpi"),
    "SRR5660044": (os.path.join(sample_data_dir, "sample_SRR5660044_1.fastq"), os.path.join(sample_data_dir, "sample_SRR5660044_2.fastq"), "Donor 3", "2dpi"),
    "SRR5660045": (os.path.join(sample_data_dir, "sample_SRR5660045_1.fastq"), os.path.join(sample_data_dir, "sample_SRR5660045_2.fastq"), "Donor 3", "6dpi")
}

# Create output directory for Bowtie2 results
output_dir = os.path.join(pipeline_dir, "Bowtie2_output")
os.makedirs(output_dir, exist_ok=True)

# Map reads to the HCMV genome using Bowtie2 with optimizations
# Loop through the samples
for sample, (fastq1, fastq2, donor, condition) in samples.items():
    if not os.path.exists(fastq1) or not os.path.exists(fastq2):
        print(f"Error: One or both FASTQ files for sample {sample} do not exist.")
        continue  # Skip this sample and move to the next one

    # Count total reads before Bowtie2 mapping
    total_reads_command1 = f"wc -l {fastq1}"
    total_reads1_output = subprocess.check_output(total_reads_command1, shell=True).strip().decode('utf-8')  # Decode bytes to string
    total_reads1 = int(total_reads1_output.split()[0]) // 4  # Each read is 4 lines in the fastq file

    total_reads_command2 = f"wc -l {fastq2}"
    total_reads2_output = subprocess.check_output(total_reads_command2, shell=True).strip().decode('utf-8')
    total_reads2 = int(total_reads2_output.split()[0]) // 4
    total_reads = total_reads1 + total_reads2

    # Run Bowtie2
    bowtie2_command = (
        f"bowtie2 --very-fast-local -x {index_prefix} -1 {fastq1} -2 {fastq2} "
        f"-S {os.path.join(output_dir, f'{sample}.sam')} --quiet"
    )
    try:
        subprocess.run(bowtie2_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Bowtie2 for sample {sample}: {e}")
        continue  # Skip to the next sample if Bowtie2 fails

    # Convert SAM to BAM and filter mapped reads
    sam_file = os.path.join(output_dir, f"{sample}.sam")
    bam_file = os.path.join(output_dir, f"{sample}.bam")
    sorted_bam_file = os.path.join(output_dir, f"{sample}_sorted.bam")
    mapped_bam_file = os.path.join(output_dir, f"{sample}_mapped.bam")

    if not os.path.exists(sam_file):
        print(f"Error: SAM file for sample {sample} does not exist.")
        continue  # Skip this sample and move to the next one

    # Convert SAM to BAM
    samtools_view_command = f"samtools view -bS {sam_file} > {bam_file}"
    try:
        subprocess.run(samtools_view_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error converting SAM to BAM for sample {sample}: {e}")
        continue  # Skip to the next sample if SAM to BAM conversion fails

    # Sort BAM file
    samtools_sort_command = f"samtools sort {bam_file} -o {sorted_bam_file}"
    subprocess.run(samtools_sort_command, shell=True, check=True)

    # Index BAM file
    samtools_index_command = f"samtools index {sorted_bam_file}"
    subprocess.run(samtools_index_command, shell=True, check=True)

    # Filter mapped reads
    samtools_filter_command = f"samtools view -b -F 4 {sorted_bam_file} > {mapped_bam_file}"
    subprocess.run(samtools_filter_command, shell=True, check=True)

    # Count reads after filtering
    mapped_reads_command = f"samtools view -c {mapped_bam_file}"
    mapped_reads = int(subprocess.check_output(mapped_reads_command, shell=True).strip())

    # Write Bowtie2 log
    write_bowtie2_log(log_file_path, donor, condition, total_reads, mapped_reads)

print("Bowtie2 mapping and filtering complete.")