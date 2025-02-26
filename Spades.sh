#!/bin/bash

# Log file for SPAdes commands
log_file="/home/2025/mbrady9/PipelineProject_Michaela_Brady/PipelineProject.log"

# Log the start of the process for Donor 1 (2dpi)
echo "Running SPAdes for Donor 1 (2dpi)" >> $log_file
echo "spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660030_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660030_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor1_2dpi --careful --threads 4" >> $log_file

# Run SPAdes for Donor 1 (2dpi)
spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660030_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660030_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor1_2dpi --careful --threads 4

# Log the start of the process for Donor 1 (6dpi)
echo "Running SPAdes for Donor 1 (6dpi)" >> $log_file
echo "spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660033_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660033_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor1_6dpi --careful --threads 4" >> $log_file

# Run SPAdes for Donor 1 (6dpi)
spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660033_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660033_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor1_6dpi --careful --threads 4

# Log the start of the process for Donor 3 (2dpi)
echo "Running SPAdes for Donor 3 (2dpi)" >> $log_file
echo "spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660045_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660045_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor3_2dpi --careful --threads 4" >> $log_file

# Run SPAdes for Donor 3 (2dpi)
spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660045_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660045_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor3_2dpi --careful --threads 4

# Log the start of the process for Donor 3 (6dpi)
echo "Running SPAdes for Donor 3 (6dpi)" >> $log_file
echo "spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660044_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660044_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor3_6dpi --careful --threads 4" >> $log_file

# Run SPAdes for Donor 3 (6dpi)
spades.py -1 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660044_1.fastq -2 /home/2025/mbrady9/PipelineProject_Michaela_Brady/sample_data/sample_SRR5660044_2.fastq -k 77 -o /home/2025/mbrady9/PipelineProject_Michaela_Brady/Donor_spades/Donor3_6dpi --careful --threads 4

# Final log message
echo "SPAdes assembly for all donors and timepoints completed." >> $log_file
