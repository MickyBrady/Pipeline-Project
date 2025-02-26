import os
from Bio import SeqIO

# Path to the SPAdes assemblies
spades_dir = "/home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades"
log_file_path = "/home/2025/mbrady9/PipelineProject_Michaela_Brady/PipelineProject.log"

# Function to find the longest contig in each SPAdes assembly
def get_longest_contig(assembly_file):
    longest_contig = None
    max_length = 0
    
    for record in SeqIO.parse(assembly_file, "fasta"):
        contig_length = len(record.seq)
        if contig_length > max_length:
            longest_contig = record
            max_length = contig_length
    
    return longest_contig

# Loop through the subdirectories and retrieve the longest contig for each assembly
longest_contigs = {}
for subdir in os.listdir(spades_dir):
    subdir_path = os.path.join(spades_dir, subdir)
    if os.path.isdir(subdir_path):
        # Loop through each FASTA file in the subdirectory
        for filename in os.listdir(subdir_path):
            if filename.endswith(".fasta"):  # Assuming the assembly files are in FASTA format
                assembly_file = os.path.join(subdir_path, filename)
                longest_contig = get_longest_contig(assembly_file)
                if longest_contig:
                    longest_contigs[f"{subdir}_{filename}"] = longest_contig

# Save the longest contigs to files and keep a record
output_dir = "/home/2025/mbrady9/PipelineProject_Michaela_Brady/Longest_Contigs"
os.makedirs(output_dir, exist_ok=True)

for sample, contig in longest_contigs.items():
    output_file = os.path.join(output_dir, f"{sample}_longest_contig.fasta")
    with open(output_file, "w") as output_handle:
        SeqIO.write(contig, output_handle, "fasta")
    print(f"Sample: {sample}")
    print(f"Longest Contig: {contig.id}")
    print(f"Length: {len(contig.seq)}")
    print(f"Sequence: {contig.seq[:100]}...")  # Print the first 100 bases
    print(f"Saved to: {output_file}")

# Retrieve Betaherpesvirinae complete genomes from NCBI
os.system('esearch -db nucleotide -query "Betaherpesvirinae[Organism] AND complete genome[Title] AND srcdb_refseq[PROP]" | efetch -format fasta > betaherpes_refseq.fasta')

# Example usage for BLAST+ query
# Create a local BLAST database using the provided betaherpes_refseq.fasta file
def create_blast_db(fasta_file, db_name):
    os.system(f"makeblastdb -in {fasta_file} -dbtype nucl -out {db_name}")

# Run BLAST+ and keep the best alignment
def run_blast(query_file, db_name, output_file):
    os.system(f"blastn -query {query_file} -db {db_name} -out {output_file} -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -max_target_seqs 10")

# Create BLAST database
create_blast_db("betaherpes_refseq.fasta", "betaherpes_db")

# Run BLAST+ for each longest contig and write results to log file
with open(log_file_path, "a") as log_file:
    for sample, contig in longest_contigs.items():
        query_file = os.path.join(output_dir, f"{sample}_longest_contig.fasta")
        output_file = os.path.join(output_dir, f"{sample}_blast_results.txt")
        run_blast(query_file, "betaherpes_db", output_file)
        
        donor = "Donor1" if "Donor1" in sample else "Donor3"
        log_file.write(f"{donor}:\n")
        log_file.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
        with open(output_file, "r") as blast_results:
            for line in blast_results:
                log_file.write(line)

print("BLAST+ queries completed. Results are saved in the log file.")