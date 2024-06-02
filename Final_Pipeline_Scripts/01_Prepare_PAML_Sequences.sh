#!/bin/bash

# This script requires the following variables be 'exported' into the
# job script envrionment:
#	_PIPE_SCRIPTS_FROM_GITHUB
#	_PIPE_SCRATCH_DIR
#	_PIPE_FINAL_OUTPUT_DIR
#	_PIPE_COHORT_MEMBERS
#	_PIPE_RUN_NICKNAME
#	_PIPE_ALL_DATA

# Look for the checkpoint. Exit with success if we find it, exit without error.
if [ -f "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/01_Prepare_PAML_Sequences.done" ]
then
    exit 0
fi
# Set bash "strict mode"
set -euo pipefail

# Load conda environment
module load anaconda
conda activate CYPevol

# Define paths to the Python scripts for translating, backtranslating, and selecting
# representative orthologues. These are distributed via GitHub repository, so the
# user will have a local copy when they clone the repo.
TRANSLATE_SCRIPT="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/Translate_Orthogroup_Sequences.py"
BACKTRANSLATE_SCRIPT="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/Backtranslate_AA_aligned.py"
REP_ORTHOLOGUE_SCRIPT="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/Choose_Representative_Orthologs.py"

# Define path to directory of orthogroup sequences from Orthofinder. This is the directory from the
# Orthofinder script where the orthogroups with human CYPs of interest were copied into.
TARGET_CYP_DIR_FULLPATH="${_PIPE_FINAL_OUTPUT_DIR}/Step_00_Orthofinder_TargetOGs"

# Define path to the directory with CDS FASTA files that were used as Orthofinder inputs
# We will identify the full path to the first entry of the list of species, and then use its
# directory name (dirname) to get the path to where the CDS FASTAs are deposited.
ORTHOFINDER_CDS_INPUT_DIR="$(dirname "${_PIPE_COHORT_MEMBERS[0]}")"

# Define path to the Orthogroups.tsv from Orthofinder
ORTHOGROUPS_TSV="${_PIPE_ALL_DATA}/Orthogroups_${_PIPE_RUN_NICKNAME}.tsv"

# Define the path to the CYP-Protein ID-Orthogroup ID CSV file
TARGET_CYP_CSV="${_PIPE_ALL_DATA}/CYPnames_Trans_Prot_with_OGs_${_PIPE_RUN_NICKNAME}.csv"

# Define paths to hold the intermediate files.
# A path for orthogroup sequences that have been translated to amino acids
OG_AA_DIR="${_PIPE_SCRATCH_DIR}/${_PIPE_RUN_NICKNAME}_Target_OG_AA_Seqs"
# A path for the aligned amino acid orthogroup sequences
OG_AA_ALIGN_DIR="${_PIPE_SCRATCH_DIR}/${_PIPE_RUN_NICKNAME}_Target_OG_AA_Align_Seqs"
# A path for the backtranslated aligned orthogroup sequences
OG_BACKTRANSLATED_ALIGN_DIR="${_PIPE_SCRATCH_DIR}/${_PIPE_RUN_NICKNAME}_Target_OG_Backtranslated_Align_Seqs"
# A path to the backtranslated, aligned, one representative sequence per species orthogroup sequences. These will be PAML inputs
PAML_SEQ_INPUT_DIR=""${_PIPE_FINAL_OUTPUT_DIR}/Step_01_PAML_Seq_Inputs""


# Make the output directories if they do not yet exist on disk
mkdir -p "${OG_AA_DIR}" "${OG_AA_ALIGN_DIR}" "${OG_BACKTRANSLATED_ALIGN_DIR}" "${PAML_SEQ_INPUT_DIR}"

# Run the translate-align-backgranslate-representative orthologue steps in a for loop
for OG_FASTA in $(find "${TARGET_CYP_DIR_FULLPATH}" -mindepth 1 -maxdepth 1 -type f -name '*.fa')
do
	# Extract the OG ID from the filename
	OG_ID=$(basename "${OG_FASTA}" | sed -e 's/.fa//g')
	#	First, translate the orthogroup sequence to amino acid
	python "${TRANSLATE_SCRIPT}" "${OG_FASTA}" > "${OG_AA_DIR}/${OG_ID}_AA.fa"

	# 	Second, align the amino acids with MAFFT-linsi
	mafft-linsi --maxiterate 1000 "${OG_AA_DIR}/${OG_ID}_AA.fa" > "${OG_AA_ALIGN_DIR}/${OG_ID}_AA_Aligned.fa"

	#	Third, backtranslate the aligned amino acids to nucleotide
	python "${BACKTRANSLATE_SCRIPT}" "${ORTHOFINDER_CDS_INPUT_DIR}" "${OG_AA_ALIGN_DIR}/${OG_ID}_AA_Aligned.fa" "${ORTHOGROUPS_TSV}" > "${OG_BACKTRANSLATED_ALIGN_DIR}/${OG_ID}_Backtranslated.fa"

	#	Fourth, select one representative orthologue from each species, fix names for PAML, and replace gap (-) with question mark (?)
	#		Note that for humans, the representative will be the NCBI sequence that is officially associated with that CYP gene
	python "${REP_ORTHOLOGUE_SCRIPT}" "${OG_BACKTRANSLATED_ALIGN_DIR}/${OG_ID}_Backtranslated.fa" "${TARGET_CYP_CSV}" > "${PAML_SEQ_INPUT_DIR}/${OG_ID}_RepOrthologues.fa"
done

# Put a checkpoint file when we finish
touch "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/01_Prepare_PAML_Sequences.done"
