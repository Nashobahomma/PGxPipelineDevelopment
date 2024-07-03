#!/bin/bash
# Authors:
#   Amber Nashoba
#   Tom Kono
# This script requires the following variables be 'exported' into the
# job script envrionment:
#	_PIPE_FINAL_OUTPUT_DIR
# 	_PIPE_SCRATCH_DIR
#	_PIPE_RUN_NICKNAME
#	_PIPE_ALL_DATA


# Look for the checkpoint. Exit with success if we find it, exit without error.
if [ -f "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/04_Run_PAML_Site_Models.done" ]
then
    exit 0
fi
# Set bash "strict mode"
set -euo pipefail

# Load conda environment
module load anaconda
conda activate CYPevol

# We will put a "codeml working directory" in /scratch because codeml writes ~10 files during the course of its
# analysis (in addition to the main .out file). These will be in scratch because we will potentially be generating
# hundreds of these little files, which is not good for /projects.
CODEML_WORKING_DIRECTORY="${_PIPE_SCRATCH_DIR}/${_PIPE_RUN_NICKNAME}_codeml_work_dir"
# Define a path to where to copyt the "codeml accessory" files. These are not the main .out files, but are the
# rub, rst, etc.  that are made during PAML runs. This directory is the final destination of the zipped
# accessory files from each run.
CODEML_ACCESSORY_OUT="${_PIPE_FINAL_OUTPUT_DIR}/Step_04Acc_PAML_Accessory_Files"
# Define a path to the input directory of codeml control files
CONTROL_FILE_DIRECTORY="${_PIPE_FINAL_OUTPUT_DIR}/Step_03_PAML_Control_Files"

# Define path to the Python script that parses PAML output to produce a table that is
# easy to analyze with R
PARSE_PAML_OUT_PY="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/Parse_PAML_Outputs.py"

# Define a path to the R script that performs PAML model comparisons
PAML_MODEL_COMP_R="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/PAML_Site_Model_Compare_and_Summarize.R"
# Define output for the "full" PAML model comparison and "digest" PAML model comparison CSVs
FULL_PAML_MODEL_COMPARISON_CSV="${_PIPE_ALL_DATA}/Parsed_PAML_Output_table_With_Model_Comparisons_${_PIPE_RUN_NICKNAME}.csv"
DIGEST_PAML_MODEL_COMPARISON_CSV="${_PIPE_ALL_DATA}/Parsed_PAML_Output_table_Omegas_and_PValues_${_PIPE_RUN_NICKNAME}.csv"


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

# Parse the directory of PAML outputs to produce a table with dN/dS (omegas) and maximum likelihood values
# for ease of analysis in R
python "${PARSE_PAML_OUT_PY}" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_04_PAML_Runs" \
	"${_PIPE_ALL_DATA}/CYPnames_Trans_Prot_with_OGs_${_PIPE_RUN_NICKNAME}.csv" \
	> "${_PIPE_ALL_DATA}/Parsed_PAML_Output_table_${_PIPE_RUN_NICKNAME}.csv"

# Call the R script to run the model comparisons and produce the CSVs with model comparison results
Rscript "${PAML_MODEL_COMP_R}" "${_PIPE_ALL_DATA}/Parsed_PAML_Output_table_${_PIPE_RUN_NICKNAME}.csv" "${FULL_PAML_MODEL_COMPARISON_CSV}" "${DIGEST_PAML_MODEL_COMPARISON_CSV}"

# Make a checkpoint file
touch "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/04_Run_PAML_Site_Models.done"
