#!/bin/bash
# The "lemma" -- the outer layer of a seed coat in a grass flower. This is the
# script that users should edit to customize their individual pipeline runs.
# This is the configuration script for the pipeline.

# Use Bash "strict mode"
#       -e: Exit on first error, rather than try to run the whole script
#       -u: Make undefined variables errors, rather than using empty strings as default value
#       -o pipefail: when one command in a pipe fails, treat all commands as having failed
set -e
set -u
set -o pipefail
IFS=$'\n\t'

# Unset the configuration variables before we set them in this script. This is
# to prevent any "cross-run" memory. This is similar to the "rm(list=ls())"
# statement in R.
unset _PIPE_SCRIPTS_FROM_GITHUB _PIPE_ANALYSIS_START_DATE _PIPE_PARTITION \
	_PIPE_MEM_PER_CPU _PIPE_WALLTIME _PIPE_CPUS_PER_TASK _PIPE_NTASKS \
	_PIPE_NNODES _PIPE_SLURM_ACCOUNT _PIPE_EMAIL_TYPES _PIPE_SCRATCH_DIR \
	_PIPE_ALL_DATA _PIPE_ALL_CDS _PIPE_CYP_NAME_PROTEIN_ID \
	_PIPE_COHORT_MEMBERS _PIPE_RUN_NICKNAME

# Define the path to the user-specific copy of the GitHub repository. Each
# user should have their own version of the pipeline scripts. This is the path
# to where the the user has cloned the GitHub repository onto Alpine (or whichever cluster).
export _PIPE_SCRIPTS_FROM_GITHUB="/projects/anashoba@xsede.org/PipelineFiles/PGxPipelineDevelopment"

# Define the current run date as a variable that we can reference later in the script, so that if a run
# crosses the midnight time, then the dates will be consistent. This is for making sure that all ouput
# files that are named with a date stamp have consistent dates.
export _PIPE_ANALYSIS_START_DATE=$(date +'%F %T')

# Define the scheduler resource request parameters for each sbatch command. These
# should be tailored to the cluster, dataset size, and partitions available.
# These parameters should be defined for a single job because each job will be
# send with these resource requests.
#	Note that these terms are the "official" terms that the Slurm documentation
# 	uses for job resource requests
# 		Partition: The queue to which jobs should be submitted
#		Walltime: "real time"  - time as human would experience it (as opposed to CPU-hours)
#		N. tasks: The number of "sub-jobs" that each job script has. Should just stay at 1.
#		N. nodes: Number of physically separate computers in the cluster that should allocated. Keep this at 1.
#		Slurm account: the group for tracking resource usage. Set this to the Claw Lab group name on the cluster.
export _PIPE_PARTITION="amilan"
export _PIPE_MEM_PER_CPU="4gb"
export _PIPE_WALLTIME="12:00:00"
export _PIPE_CPUS_PER_TASK="16"
export _PIPE_NTASKS="1"
export _PIPE_NNODES="1"
export _PIPE_SLURM_ACCOUNT=""
export _PIPE_EMAIL_TYPES="BEGIN,END,FAIL"

# Define a path to the scratch directory on Alpine
export _PIPE_SCRATCH_DIR="/scratch/alpine/${USER}"
# Define a path to where the input data for the pipeline will be:
#	- Directory with CDS FASTA files
#	- CYP name -> Protein ID CSV
# This directory should be kept in scratch space.
export _PIPE_ALL_DATA="/projects/anashoba@xsede.org/PipelineFiles"
# Path to directory containing all cds fasta files
export _PIPE_ALL_CDS="${_PIPE_ALL_DATA}/NCBI_CDS_Files"
# Path to the table connecting CYP450 gene names with their associated the human protein IDs
export _PIPE_CYP_NAME_PROTEIN_ID="${_PIPE_ALL_DATA}/CYP57genesHumanProteinIDs.csv"

#Cohort Designation
# In this directory RENAME IT FOR POSTERITY List out the full genus_species name for each species' CDS 
# to be included in this run of the pipeline.
# For _PIPE_COHORT_MEMBERS, you must specify the full paths to the names of the
# species that you want to analyze. 
export _PIPE_COHORT_MEMBERS=(
	"${_PIPE_ALL_CDS}/Callithrix_jacchus.fasta"
	"${_PIPE_ALL_CDS}/Carlito_syrichta.fasta"
	"${_PIPE_ALL_CDS}/Gorilla_gorilla.fasta"
	"${_PIPE_ALL_CDS}/Homo_sapiens.fasta"
	"${_PIPE_ALL_CDS}/Lemur_catta.fasta"
	"${_PIPE_ALL_CDS}/Macaca_mulatta.fasta"
	"${_PIPE_ALL_CDS}/Microcebus_murinus.fasta"
	"${_PIPE_ALL_CDS}/Nomascus_leucogenys.fasta"
	"${_PIPE_ALL_CDS}/Nycticebus_coucang.fasta"
	"${_PIPE_ALL_CDS}/Otolemur_garnettii.fasta"
	"${_PIPE_ALL_CDS}/Pan_troglodytes.fasta"
	"${_PIPE_ALL_CDS}/Papio_anubis.fasta"
	"${_PIPE_ALL_CDS}/Piliocolobus_tephrosceles.fasta"
	"${_PIPE_ALL_CDS}/Pongo_abelii.fasta"
	"${_PIPE_ALL_CDS}/Rhinopithecus_roxellana.fasta"
	"${_PIPE_ALL_CDS}/Sapajus_apella.fasta")
export _PIPE_COHORT_MEMBER_NUMBER="${#_PIPE_COHORT_MEMBERS[@]}"

# Define a nickname for the analysis. This will be used to name the ouput
# directory for the pipeline.
export _PIPE_RUN_NICKNAME="${_PIPE_COHORT_MEMBER_NUMBER}_${_PIPE_ANALYSIS_START_DATE}"

# Version identifier for the pipeline.
export _PIPE_VERSION="0.0.0"
# Information about who is running the pipeline and when. Do not edit these.
export _PIPE_RUN_USER="$(id -un)"
export _PIPE_RUN_TIME="$(date '+%Y-%m-%d_%T')"

#copy directory containing all data and scripts into the scratch directory
#______NOTE: NEED TO MAKE THIS DIRECTORY OR REPLACE ForAlpine WITH A NEW DESCRIPTIVE NAME like PipelineFiles
#cp -r /projects/anashoba@xsede.org/PIPELINEFILES /scratch/alpine/$USER


#!the og.tsv and gene-protein-og table will be added in the orthofinder batch script


# MAYBE MAKE ALL THE STEP_XX directories here?


#Set the directory for all scripts
# ALL_SCRIPTS="/scratch/alpine/anashoba@xsede.org/SNAPPY_SCRIPTS"

# # NOTE NOTE NOTE is there a trick to batching several scripts that are themselves batch running code? 
# #Run bash scripts in this order:
# ORTHOFINDER_SH="{ALL_SCRIPTS}/Orthofinder.sh"
# TAB_SH="{ALL_SCRIPTS}/Translate_align_backtranslate.sh"
# GENE_TREES="{ALL_SCRIPTS}/Make_Gene_Trees.sh"
# CONTROL_FILES="{ALL_SCRIPTS}/Make_CtlFiles.sh"
# RUN_PAML="{ALL_SCRIPTS}/Run_Codeml_SiteModels.sh"




# THE REAL GUTS GO HERE




# #When done with the batch scripts:
# ... PARSER python script here 
# MAKE_PAML_SITEMODEL_RESULTS_TABLE="{ALL_SCRIPTS}/PAML_Output_Parser.py"
# python ""