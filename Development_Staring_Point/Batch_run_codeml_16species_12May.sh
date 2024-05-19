#!/bin/bash

# JOB RESOURCES
###############

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=4gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=16:00:00 
#SBATCH --mail-user=amber.nashoba@cuanschutz.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e PAML_12May_%j.err
#SBATCH -o PAML_12May_%j.out


# JOB COMMANDS
module load anaconda
conda activate CYPevol


cp -r /projects/anashoba@xsede.org/ForAlpine/Step_03_PAML_Seq_Inputs /scratch/alpine/$USER
cp -r /projects/anashoba@xsede.org/ForAlpine/Step_04_PAML_Gene_Trees /scratch/alpine/$USER

# Reminder: this is the output folder for PAML, as defined in the control files.
# /projects/anashoba@xsede.org/ForAlpine/Step_06_PAML_Runs

# Run 'codeml' in batch for all target CYP-containing orthogroups. In addition
# to running 'codeml,' we will include cleanup steps so that subsequent
# 'codeml' runs do not clobber each other (i.e. each time codeml is run for a gene,
## the subsequent output is not overwriting the output from a gene run previously.)
# Usage:
#	sbatch /path/to/Batch_run_codeml.sh
# Updated 2024-04-28
# Amber Nashoba and Tom Kono


# Path to the codeml program
#NEEDS ADJUSTING: 
#CODEML="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PAML/paml4.8/bin/codeml"
# Path to the main codeml output directory
#CODEML_WORKING_DIRECTORY="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_06_PAML_Runs"
# We will put a "codeml working directory" in /scratch because codeml writes ~10 files during the course of its
# analysis (in addition to the main .out file). These will be in scratch because we will potentially be generating
# hundreds of these little files, which is not good for /projects.
CODEML_WORKING_DIRECTORY="/scratch/alpine/$USER/codeml_working_dir"
# Define a path to where to copyt the "codeml accessory" files. These are not the main .out files, but are the
# rub, rst, etc.  that are made during PAML runs. This directory is the final destination of the zipped
# accessory files from each run.
CODEML_ACCESSORY_OUT="/projects/anashoba@xsede.org/ForAlpine/Step_06_PAML_Accessories"
mkdir -p "${CODEML_ACCESSORY_OUT}"
# Define a path to the input directory of codeml control files
#CONTROL_FILE_DIRECTORY="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_05_PAML_Ctl_Files"
CONTROL_FILE_DIRECTORY="/projects/anashoba@xsede.org/ForAlpine/Step_05_PAML_Ctl_Files/"
# Now, run codeml and the cleanup steps.
for ctl_file in $(find "${CONTROL_FILE_DIRECTORY}" -mindepth 1 -maxdepth 1 -type f -name '*.txt' | sort -V)
do
	# Extract the orthogroup ID and codeml model(s) name from the control file name
	og_id_and_codeml_model=$(basename "${ctl_file}" | cut -f 1-2 -d '_')
	# Make a directory in which we will run codeml and let it make its "mess"
	mkdir -p "${CODEML_WORKING_DIRECTORY}/${og_id_and_codeml_model}"
	cd "${CODEML_WORKING_DIRECTORY}/${og_id_and_codeml_model}"
	# Run codeml, redirect the text printed to the terminal to files on disk
	#	Not sure if these will be useful to save, but we will do it just in case it is.
	# codeml seems to have a problem with long path names to control files, so we will
	# copy the control file into the codeml "working directory" and run it from there. It is
	# annoying to have to do this, but codeml is a "special boy" and needs a lot of handholding.
	# Then, to avoid having duplicate files on disk, we will remove the copied control file.
	cp "${ctl_file}" ./
	# Write some progress messages to the terminal so that we can keep track of progress
	echo "$(date '+%F %T') Running codeml on ${og_id_and_codeml_model} ..."  > /dev/stderr
	codeml "./$(basename ${ctl_file})" > "${og_id_and_codeml_model}.stdout.txt" 2> "${og_id_and_codeml_model}.stderr.txt"
	echo "$(date '+%F %T') Done running ${og_id_and_codeml_model}"  > /dev/stderr
	rm "./$(basename ${ctl_file})"
	# Cd "up" one level to zip the codeml aftermath (i.e. the extra files)
	cd ../
	zip -r "${og_id_and_codeml_model}.zip" "./${og_id_and_codeml_model}/"
	# Delete the working direcotry
	rm -rf "./${og_id_and_codeml_model}/"
	# Copy the .zip into the codeml accessory folder in /projects.
	cp "${og_id_and_codeml_model}.zip" "${CODEML_ACCESSORY_OUT}/"
done
