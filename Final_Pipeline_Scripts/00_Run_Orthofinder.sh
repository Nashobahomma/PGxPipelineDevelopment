#!/bin/bash


# This script requires the following variables be 'exported' into the
# job script envrionment:
#	_PIPE_SCRIPTS_FROM_GITHUB
#	_PIPE_SCRATCH_DIR
#	_PIPE_FINAL_OUTPUT_DIR
#	_PIPE_COHORT_MEMBERS
#	_PIPE_RUN_NICKNAME
#	_PIPE_CYP_NAME_PROTEIN_ID
#	_PIPE_ALL_DATA

# Look for the checkpoint. Exit with success if we find it, exit without error.
if [ -f "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/00_Run_Orthofinder.done" ]
then
    exit 0
fi

# JOB COMMANDS
module load anaconda
conda activate CYPevol

# Use Bash "strict mode"
#       -e: Exit on first error, rather than try to run the whole script
#       -u: Make undefined variables errors, rather than using empty strings as default value
#       -o pipefail: when one command in a pipe fails, treat all commands as having failed
set -euo pipefail

# Path to ParseOrthogroupsFromResults.py, to isolate orthogroups of interest based on having a CYP.
PARSE_ORTHOGROUPS_PY="${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/ParseOrthogroupsFromResults.py"

# Use the mktemp comand to return a directory path in the system temp space
cd "${_PIPE_SCRATCH_DIR}"
TEMP_DIR_NAME=$(mktemp -d)
cd "${TEMP_DIR_NAME}"

# Read the _PIPE_COHORT_MEMBERS array varaible to identify the species that
# the user wants to analyze. Iterate over the elements of _PIPE_COHORT_MEMBERS
# and link each file into a temp input directory
mkdir -p ./orthofinder_cds_in
for species_fasta in "${_PIPE_COHORT_MEMBERS[@]}"
do
	ln -s "${species_fasta}" "${TEMP_DIR_NAME}/orthofinder_cds_in/$(basename "${species_fasta}")"
done

# Run Orthofinder
#	Define a output folder name for Orthofinder based on the date that the run
#	is performed.
HPC_ORTHOFINDER_OUT="${_PIPE_RUN_NICKNAME}"
#	Remove the output folder if already exists, to prevent "cross-contamination" of
#	Orthofinder analyses.
rm -rf "./${HPC_ORTHOFINDER_OUT}"
echo "$(date +'%F %T'): Starting Orthofinder"
orthofinder -f "${TEMP_DIR_NAME}/orthofinder_cds_in" -t "${SLURM_CPUS_PER_TASK}" -a "${SLURM_CPUS_PER_TASK}" -d -M "msa" -os -S "blast_nucl" -A "mafft" -T "raxml" -z -o "${HPC_ORTHOFINDER_OUT}" 
echo "$(date +'%F %T'): Finished Orthofinder"


# Use a find command to get the full path to the Orthogroups.tsv file in the Orthofinder results directory
# This is because Orthofinder names its output based on the date it is run, which will necessarily change each run.
OG_TSV=$(find "${HPC_ORTHOFINDER_OUT}" -type f -name 'Orthogroups.tsv')
# Use a find command to get the path to the Orthogroup_Sequences directory in the Orthofinder results directory.
OG_SEQS_DIR=$(find "${HPC_ORTHOFINDER_OUT}" -type d -name 'Orthogroup_Sequences')
# Make the full path to the destination directory for target CYPs. Use a find command to get the name of the
# date-specific results directory that is made by Orthofinder.
ORTHOFINDER_DATE_DIR=$(find "${HPC_ORTHOFINDER_OUT}" -mindepth 1 -maxdepth 1 -type d -name 'Results_*')
TARGET_CYP_DIR_FULLPATH="${_PIPE_FINAL_OUTPUT_DIR}/Step_00_Orthofinder_TargetOGs"
python "${PARSE_ORTHOGROUPS_PY}" "${_PIPE_CYP_NAME_PROTEIN_ID}" "${OG_TSV}" "${OG_SEQS_DIR}" "${TARGET_CYP_DIR_FULLPATH}" > "${_PIPE_ALL_DATA}/CYPnames_Trans_Prot_with_OGs_${_PIPE_RUN_NICKNAME}.csv"
cp "${OG_TSV}" "${_PIPE_ALL_DATA}/Orthogroups_${_PIPE_RUN_NICKNAME}.tsv"

# Remove the temp directory that we used for analysis. Not strictly necessary,
# but nice for multiuser systems.
cd
rm -rf "${TEMP_DIR_NAME}"

# If we make to this step, then orthofinder and the parsing of orthofinder
# outputs have finished without error. Make a checkpoint file so that we know
# that this step is complete.
touch "${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints/00_Run_Orthofinder.done"
