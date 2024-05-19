#!/bin/bash
# 2023-02-26 # January 2024 -> 16 Species_57genes version and adjusted pathways
# AN and TK
# Shell script to run the translate-align-backtranslate-RepOrthologSelection steps for
# a batch of orthogroups with target CYP genes of interest. Depends on MAFFT being
# installed.


# Define paths to the Python scripts for translating, backtranslating, and selecting
# representative orthologues.
TRANSLATE_SCRIPT="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/Scripts_16species_57Genes_versions_January2024_OFruns/Translate_Orthogroup_Sequences.py"
BACKTRANSLATE_SCRIPT="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/Scripts_16species_57Genes_versions_January2024_OFruns/Backtranslate_AA_aligned.py"
REP_ORTHOLOGUE_SCRIPT="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/Scripts_16species_57Genes_versions_January2024_OFruns/OrthologSelection_16species.py"

# Define path to directory of orthogroup sequences from Orthofinder
#	This will have to be manually downloaded from GCS bucket after running Orthofinder. This directory structure /naming does not account for different run dates of OF
OGS_INPUT_DIR="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/TargetCYP_Sequences/TargetCYP_Sequences/TargetCYP_Sequences"


# Define path to the directory with CDS FASTA files that were used as Orthofinder inputs
#	This is the same directory for the "prepping to run Orthofinder" step in the pipeline. CHECK THIS
ORTHOFINDER_CDS_INPPUT_DIR="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/CDS_FASTA_16species"

# Define path to the Orthogroups.tsv from Orthofinder
#	This will have to be manually downloaded from GCS after Orthofinder is finished on HPC.
## MAKE NOTE OF THIS IN THE LAB ARCHIVE AMBER
ORTHOGROUPS_TSV="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/2024-04-07_Orthogroups_AllOGs_16Species_57Genes.tsv"

# Define the path to the CYP-Protein ID-Orthogroup ID CSV file
#	This is downloaded from GCS after Orthofinder is finished on HPC. PROBLEM: THIS TABLE HAS ERRORS WHEN MULTIPLE CYPGENE-HUMAN PROTEINS ARE FOUND IN THE SAME ORTHOGROUP
#TARGET_CYP_CSV="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/TargetCYP_Sequences/CYPnames_Trans_Prot_with_OGs_16species_57Genes.csv"
#temptry
##Try this again - but what happens if non-unique OG-human protein name pair? 
TARGET_CYP_CSV="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/CYPnames_Trans_Prot_with_OGs_16Species_57Genes.csv"

# Define paths to hold the intermediate files.
# A path for orthogroup sequences that have been translated to amino acids
OG_AA_DIR="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_00_Target_OG_AA_Seqs"
# A path for the aligned amino acid orthogroup sequences
OG_AA_ALIGN_DIR="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_01_Target_OG_AA_Align_Seqs"
# A path for the backtranslated aligned orthogroup sequences
OG_BACKTRANSLATED_ALIGN_DIR="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_02_Target_OG_Backtranslated_Align_Seqs"
# A path to the backtranslated, aligned, one representative sequence per species orthogroup sequences. These will be PAML inputs
PAML_SEQ_INPUT_DIR="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_03_PAML_Seq_Inputs"


# Make the output directories if they do not yet exist on disk
mkdir -p "${OG_AA_DIR}" "${OG_AA_ALIGN_DIR}" "${OG_BACKTRANSLATED_ALIGN_DIR}" "${PAML_SEQ_INPUT_DIR}"

# Run the translate-align-backgranslate-representative orthologue steps in a for loop
for OG_FASTA in $(find "${OGS_INPUT_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.fa')
do
	# Extract the OG ID from the filename
	OG_ID=$(basename "${OG_FASTA}" | sed -e 's/.fa//g')
	#	First, translate the orthogroup sequence to amino acid
	python "${TRANSLATE_SCRIPT}" "${OG_FASTA}" > "${OG_AA_DIR}/${OG_ID}_AA.fa"

	# 	Second, align the amino acids with MAFFT-linsi
	mafft-linsi --maxiterate 1000 "${OG_AA_DIR}/${OG_ID}_AA.fa" > "${OG_AA_ALIGN_DIR}/${OG_ID}_AA_Aligned.fa"

	#	Third, backtranslate the aligned amino acids to nucleotide
	python "${BACKTRANSLATE_SCRIPT}" "${ORTHOFINDER_CDS_INPPUT_DIR}" "${OG_AA_ALIGN_DIR}/${OG_ID}_AA_Aligned.fa" "${ORTHOGROUPS_TSV}" > "${OG_BACKTRANSLATED_ALIGN_DIR}/${OG_ID}_Backtranslated.fa"

	#	Fourth, select one representative orthologue from each species, fix names for PAML, and replace gap (-) with question mark (?)
	#		Note that for humans, the representative will be the NCBI sequence that is officially associated with that CYP gene
	python "${REP_ORTHOLOGUE_SCRIPT}" "${OG_BACKTRANSLATED_ALIGN_DIR}/${OG_ID}_Backtranslated.fa" "${TARGET_CYP_CSV}" > "${PAML_SEQ_INPUT_DIR}/${OG_ID}_RepOrthologues.fa"
done

# Leave these comments in for reference; we do not need them at present.
# First: Call the python script: `Translate_Orthogroup_Sequences.py'
# Second: run mafft
## mafft-linsi --maxiterate 1000 ~/Desktop/Dropbox/KClaw/PipelineA/13JulyOrthofinder/Unaligned_AA_Sequences/OG0006627_AA.fa > ~/Desktop/Dropbox/KClaw/PipelineA/13JulyOrthofinder/Aligned_AA_Sequences/OG0006627_AA_Aligned.fa
# Third run `Backtranslate_AA_aligned.py`

