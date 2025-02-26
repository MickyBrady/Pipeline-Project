from Bio import Entrez, SeqIO
import os

# Set up email for NCBI Entrez
Entrez.email = 'mbrady9@luc.edu'

# Set up the accession to retrieve from NCBI
hcmv_genome_id = "NC_006273.2"
handle = Entrez.efetch(db="nucleotide", id=hcmv_genome_id, rettype="gb", retmode="text")

# Parse the GenBank file
record = SeqIO.read(handle, "genbank")

# Extract CDS features
cds_sequences = []
for feature in record.features:
    if feature.type == "CDS":
        # Extract the sequence and add it to the list
        if "protein_id" in feature.qualifiers:
            cds_seq = feature.extract(record.seq)
            cds_sequences.append((feature.qualifiers["protein_id"][0], cds_seq))

# Save CDS sequences to a FASTA file
fasta_filename = "HCMV_cds.fasta"
with open(fasta_filename, "w") as f:
    for protein_id, cds_seq in cds_sequences:
        f.write(f">{protein_id}\n")
        f.write(str(cds_seq) + "\n")

# Count the number of CDS
num_cds = len(cds_sequences)

# Set directory to store log
log_directory = "/home/2025/mbrady9/PipelineProject_Michaela_Brady"
os.makedirs(log_directory, exist_ok=True)

# Change to the log directory
os.chdir(log_directory)

# Define log file path
log_file_path = os.path.join(log_directory, "PipelineProject.log")

# Check if log message already exists
log_message = f"The HCMV genome (NC_006273.2) has {num_cds} CDS."
if os.path.exists(log_file_path):
    with open(log_file_path, "r") as log_file:
        if log_message in log_file.read():
            print("Skipping redundant HCMV log entry.")
        else:
            with open(log_file_path, "a") as log_file:
                log_file.write(log_message + "\n")
            print(log_message)  # Print only if it's logged for the first time
else:
    with open(log_file_path, "w") as log_file:
        log_file.write("Pipeline Execution Log\n")
        log_file.write(log_message + "\n")
    print(log_message)
