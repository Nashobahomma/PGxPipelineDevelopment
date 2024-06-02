#!/bin/bash
# Authors:
#   Amber Nashoba
#   Tom Kono
# This script requires the following variables be 'exported' into the
# job script envrionment:
#    _PIPE_FINAL_OUTPUT_DIR

# Look for the checkpoint. Exit with success if we find it, exit without error.
if [ -f "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/02_Make_Gene_Trees.done" ]
then
    exit 0
fi
# Set bash "strict mode"
set -euo pipefail

# Load conda environment
module load anaconda
conda activate CYPevol

# Define parameters for the RAxML run - bootstrap replicates and nucletoide substitution model
BOOTSTRAP_REPLICATES="10000"
NUC_SUBSTIUTION_MODEL="GTR+G"

# Define a path to where to read PAML sequence inputs from
PAML_SEQ_INPUT_DIR="${_PIPE_FINAL_OUTPUT_DIR}/Step_01_PAML_Seq_Inputs"
# Define a path to where to write the gene trees
TREE_OUTPUT="${_PIPE_FINAL_OUTPUT_DIR}/Step_02_PAML_Gene_Trees"

for gene_alignment in $(find "${PAML_SEQ_INPUT_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.fa')
do 
    orthogroup_id="$(basename ${gene_alignment} | cut -d '_' -f 1)"
    raxml-ng \
        --all \
        --msa "${gene_alignment}" \
        --msa-format FASTA \
        --data-type DNA \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --model "${NUC_SUBSTIUTION_MODEL}" \
        --bs-trees "${BOOTSTRAP_REPLICATES}" \
        --prefix "${TREE_OUTPUT}/${orthogroup_id}"
done


# Make a checkpoint file
touch "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/02_Make_Gene_Trees.done"
