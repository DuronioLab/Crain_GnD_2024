import os
import subprocess
import shutil
import csv

# Path to sample sheet
# Replace with appropriate path to sample sheet
sample_sheet = "sample_sheet.txt"

# Path to Salmon index
# Replace with appropriate path to salmon index
salmon_index = "salmon_index"

# Path to folder containing input fastq files
# Replace with appropriate path to fastq files
input_folder = "/Path_to_Fastq/"

# Path to output folder
output_folder = "salmon_quant/"

# Library type
library_type = "A"

# Additional SLURM options
slurm_options = "--nodes=1 --cpus-per-task=4 --mem=16GB -t 60"

#Add empty list to hold sbatch commands
sbatch_commands = []

samples = []

# Read the sample sheet and extract relevant columns, skipping the first row
with open(sample_sheet, 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        sample_info = {
            'base_filename': row['useName'],
            'r1_file': row['Read1']
        }
        samples.append(sample_info)

# Loop over all fastq files in the input folder
for sample in samples:
        base_filename = sample['base_filename']
        r1_file = os.path.join(input_folder, sample['r1_file'])
        output_dir = os.path.join(output_folder, base_filename)
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Construct the Salmon command and write to bash script
        script_filename = os.path.join(output_dir, "salmon_quant_{}.sh".format(base_filename))
        salmon_command = "salmon quant -i {} -l {} -r {} -o {} --validateMappings --seqBias --useVBOpt --numBootstraps 30".format(salmon_index, library_type, r1_file, output_dir)
        with open(script_filename, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --job-name={}\n".format(base_filename))
            f.write("#SBATCH --output={}.out\n".format(base_filename))
            f.write("#SBATCH --error={}.err\n".format(base_filename))
            f.write("#SBATCH --time=6:00:00\n")
            f.write("#SBATCH --mem=16G\n")
            f.write("#SBATCH --cpus-per-task=8\n")
            f.write("echo 'Starting job {} on:'\n".format(base_filename))
            f.write("date\n")
            f.write("echo 'Running Salmon with command:'\n")
            f.write("echo '{}'\n".format(salmon_command))
            f.write("{}\n".format(salmon_command))
            f.write("echo 'Finished job {} on:'\n".format(base_filename))
            f.write("date\n")

        # Make the sbatch script executable and set file permissions
        os.chmod(script_filename, 0o755)

        # Add the sbatch command to execute the script to a list
        sbatch_command = "sbatch {}".format(script_filename)

        # Run sbatch command
        subprocess.run(sbatch_command.split())

        
