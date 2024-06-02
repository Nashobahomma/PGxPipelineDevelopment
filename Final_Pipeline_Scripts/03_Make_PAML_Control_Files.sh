#!/bin/bash
# Authors:
#   Amber Nashoba
#   Tom Kono
# This script requires the following variables be 'exported' into the
# job script envrionment:
#    _PIPE_FINAL_OUTPUT_DIR
# 	_PIPE_SCRIPTS_FROM_GITHUB

# Look for the checkpoint. Exit with success if we find it, exit without error.
if [ -f "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/03_Make_PAML_Control_Files.done" ]
then
    exit 0
fi
# Set bash "strict mode"
set -euo pipefail

# Load conda environment
module load anaconda
conda activate CYPevol

# Define path to the directory of PAML input sequence files
PAML_SEQ_INPUT_DIR="${_PIPE_FINAL_OUTPUT_DIR}/Step_01_PAML_Seq_Inputs"
# Define path to the directory of PAML input gene trees
PAML_TREE_INPUT_DIR="${_PIPE_FINAL_OUTPUT_DIR}/Step_02_PAML_Gene_Trees"
# Define path to the output directory for the control files
CONTROL_FILE_OUTPUT_DIR="${_PIPE_FINAL_OUTPUT_DIR}/Step_03_PAML_Control_Files"
# Define a path to the PAML output directory. This variable gets written into the
# control file.
PAML_OUT_DIR="${_PIPE_FINAL_OUTPUT_DIR}/Step_04_PAML_Runs"

# Define path to the Python script that actually writes a single control file
# This script is distributed with the GitHub repository, so users will have a local
# copy when they clone the repo.
WRITE_CONTROL_FILES_PY="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/Write_Site_Model_Control_File.py"


for seq_file in $(find "${PAML_SEQ_INPUT_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.fa')
	do
	orthogroup_id="$(basename ${seq_file} | cut -d '_' -f 1)"
	python3 "${WRITE_CONTROL_FILES_PY}" \
		"${seq_file}" \
		"${GENE_TREE_INPUT_DIR}/${orthogroup_id}.raxml.bestTree" \
		8a \
		"${PAML_OUT_DIR}" \
		> "${CONTROL_FILE_OUTPUT_DIR}/${orthogroup_id}_8a_Ctl_File.txt"
	python3 "${WRITE_CONTROL_FILES_PY}" \
		"${seq_file}" \
		"${GENE_TREE_INPUT_DIR}/${orthogroup_id}.raxml.bestTree" \
		01278 \
		"${PAML_OUT_DIR}" \
		> "${CONTROL_FILE_OUTPUT_DIR}/${orthogroup_id}_01278_Ctl_File.txt"
	done


# Make a checkpoint file
touch "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/03_Make_PAML_Control_Files.done"
