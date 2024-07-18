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

# Define a function to prevent RAxML errors from killing the pipeline. At least one
# RAxML error is when a FASTA file has fewer than 4 sequences.
bypass_raxml_error() {
    OG_INPUT="${1}"
    echo "RAxML caught an error on orthgroup ${OG_INPUT}" > /dev/stderr
    echo "Check the job log for the specific error message." > /dev/stderr
    exit 0
}

for gene_alignment in $(find "${PAML_SEQ_INPUT_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.fa')
do 
    orthogroup_id="$(basename ${gene_alignment} | cut -d '_' -f 1)"
    # RAxML fails on FASTA files with less than 4 sequences. We will skip these and
    # print a message to stderr indicating that we found one
    n_seqs=$(grep -c '>' ${gene_alignment})
    if [ ${n_seqs} -lt 4 ]
    then
        echo "Fewer than 4 sequences in ${orthogroup_id}; skipping RAxML for this OG." >> /dev/stderr
    else
        raxml-ng \
            --all \
            --redo \
            --msa "${gene_alignment}" \
            --msa-format FASTA \
            --data-type DNA \
            --threads "${SLURM_CPUS_PER_TASK}" \
            --model "${NUC_SUBSTIUTION_MODEL}" \
            --bs-trees "${BOOTSTRAP_REPLICATES}" \
            --prefix "${TREE_OUTPUT}/${orthogroup_id}" || bypass_raxml_error "${gene_alignment}"
        # Starting with PAML 4.10, the tree file has a required first line that contains two digits:
        #    1: number of "species" (really, sequences) in the alignment
        #    2: number of trees (this will always be 1 for our pipeline)
        # First, copy the original .raxml.bestTree output to a backup
        cp "${TREE_OUTPUT}/${orthogroup_id}.raxml.bestTree" "${TREE_OUTPUT}/${orthogroup_id}.raxml.bestTree.orig"
        # Then, generate the new PAML input tree with the required header line
        echo "${n_seqs} 1" > "${TREE_OUTPUT}/${orthogroup_id}.raxml.bestTree"
        cat "${TREE_OUTPUT}/${orthogroup_id}.raxml.bestTree.orig" >> "${TREE_OUTPUT}/${orthogroup_id}.raxml.bestTree"
    fi
done


# Make a checkpoint file
touch "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/02_Make_Gene_Trees.done"
