#!/bin/bash
# Authors:
#   Amber Nashoba
#   Tom Kono
# The "palea" --  the inner layer of a seed coat in a grass flower. This is the
# script that does all the 'sbatch' calls in the right order to execute the
# pipeline. Users should not have to edit this script.
set -e
set -u
set -o pipefail
IFS=$'\n\t'

# Unset the pipeline variables we define in this script, lest we get pollution
# of our pipeline environment.
unset _PIPE_CTLDIR _PIPE_TOT_MEM _PIPE_FINAL_OUTPUT_DIR _PIPE_EXEC_RECORD STEP_00

# source the configuration script (Lemma.sh). We will require that it exist in the same
# directory as this script.
#	This statement will return the full path to the directory in which the current script
# 	is stored. E.g., /projects/alpine/user1/dir1/script.sh -> /projects/alpine/user1/dir1
_PIPE_CTLDIR="$(dirname "${BASH_SOURCE[0]}")"
if [ -s "${_PIPE_CTLDIR}/Lemma.sh" ]
    then source "${_PIPE_CTLDIR}/Lemma.sh"
    else
        echo "Lemma.sh should exist in the same directory as the Palea.sh" > /dev/stderr
        exit 3
fi

# Define a usage function
usage() {
cat << END 1>&2
Usage:
${0}

Run the PGx CYP450 PAML analysis. These are codeml site models at present (May
2024). Users should edit the "Lemma.sh" script that is distributed with this
script ("Palea.sh") to suit their analysis and compute environment before
running this script.
END
}

# Let's run a couple resource checks here. We would need at least 16 CPUs
# and at least 48GB of memory to run.
if [ "${_PIPE_CPUS_PER_TASK}" -lt 16 ]
then
    echo "This pipeline requires at least 16 CPUs. Please set _PIPE_CPUS_PER_TASK to 16 or higher in Lemma.sh" 1>&2
    exit 1
fi

_PIPE_TOT_MEM="$(echo "$(echo "${_PIPE_MEM_PER_CPU}" | sed -E 's/gb//g') * ${_PIPE_CPUS_PER_TASK}" | bc -l)"
if [ "${_PIPE_TOT_MEM}" -lt 48 ]
then
    echo "This pipeline requires at least 48GB of RAM. Please ensure that the product of _PIPE_CPUS_PER_TASK and _PIPE_MEM_PER_CPU is at least 48GB." 1>&2
    exit 1
fi

# Check that the specified _PIPE_COHORT_MEMBERS FASTA files can each be read
for input_cds in "${_PIPE_COHORT_MEMBERS[@]}"
do 
	if [ ! -s "${input_cds}" ]
		then echo "Error: CDS file ${input_cds} could not be found. Check the filenames in the _PIPE_COHORT_MEMBERS variable in Lemma.sh" > /dev/stderr
		exit 1
	fi
done

# Make the output directory and record who is running the pipeline and how it
# is configured
_PIPE_FINAL_OUTPUT_DIR="${_PIPE_ALL_DATA}/${		}"
# Make a text file that records who, when, and where (working directory on Alpine) the pipeline was run.
_PIPE_EXEC_RECORD="${_PIPE_FINAL_OUTPUT_DIR}/PGx_Pipeline_Execution_Record.txt"
mkdir -p "${_PIPE_FINAL_OUTPUT_DIR}"
echo "Record of PGx CPY450 PAML pipeline." > "${_PIPE_EXEC_RECORD}"
echo "Pipeline version: ${_PIPE_VERSION}" >> "${_PIPE_EXEC_RECORD}"
echo "Execution time: ${_PIPE_RUN_TIME}" >> "${_PIPE_EXEC_RECORD}"
echo "Executing user: ${_PIPE_RUN_USER}" >> "${_PIPE_EXEC_RECORD}"
echo "Nickname of run: ${_PIPE_RUN_NICKNAME}" >> "${_PIPE_EXEC_RECORD}"
echo "Working directory: $(pwd -P)" >> "${_PIPE_EXEC_RECORD}"
echo "Invoking command: ${BASH_SOURCE} $(printf ' %q' "${@}")" >> "${_PIPE_EXEC_RECORD}"
echo "Input CDS FASTA directory: ${_PIPE_ALL_CDS}" >> "${_PIPE_EXEC_RECORD}"
echo "Number of CDS FASTA files: ${#_PIPE_COHORT_MEMBERS[@]}" >> "${_PIPE_EXEC_RECORD}"
echo "Final output directory: ${_PIPE_FINAL_OUTPUT_DIR}" >> "${_PIPE_EXEC_RECORD}"

# Make the directories for each step of the pipeline
mkdir -p "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Checkpoints" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_00_Orthofinder_TargetOGs" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_01_PAML_Seq_Inputs" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_02_PAML_Gene_Trees" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_03_PAML_Control_Files" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_04_PAML_Runs" \
	"${_PIPE_FINAL_OUTPUT_DIR}/Step_08_ParseOutFiles"

# Next time: get these sbatch commands finalized as we verfiy and generalize the
# job scripts.
# Submit the pipeline jobs!
#  Step 00: Run Orhofinder
STEP_00=$(sbatch \
    --parsable \
    --mail-type="${_PIPE_EMAIL_TYPES}" \
    -J "00_Run_Orthofinder" \
    -o "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/00_Run_Orthofinder.stdout" \
    -e "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/00_Run_Orthofinder.stderr" \
    -N "${_PIPE_NNODES}" \
    -n "${_PIPE_NTASKS}" \
    -c "${_PIPE_CPUS_PER_TASK}" \
    -t "${_PIPE_WALLTIME}" \
    --mem-per-cpu "${_PIPE_MEM_PER_CPU}" \
    -p "${_PIPE_PARTITION}" \
    --export="_PIPE_SCRIPTS_FROM_GITHUB=${_PIPE_SCRIPTS_FROM_GITHUB},_PIPE_CYP_NAME_PROTEIN_ID=${_PIPE_CYP_NAME_PROTEIN_ID},_PIPE_SCRATCH_DIR=${_PIPE_SCRATCH_DIR},_PIPE_RUN_NICKNAME=${_PIPE_RUN_NICKNAME},_PIPE_ALL_DATA=${_PIPE_ALL_DATA},_PIPE_FINAL_OUTPUT_DIR=${_PIPE_FINAL_OUTPUT_DIR},_PIPE_COHORT_MEMBERS=${_PIPE_COHORT_MEMBERS}" \
    "${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/00_Run_Orthofinder.sh")
echo "Step 00: Run_Orthofinder has job ID ${STEP_00}" | tee -a "${_PIPE_EXEC_RECORD}"

# Step 01: Translate, align, backtranslate, select representative orthologs
# from orthogroups with human CYPs of interest
STEP_01=$(sbatch \
    --parsable \
    --kill-on-invalid-dep=yes \
    --dependency=afterok:"${STEP_00}" \
    --mail-type="${_PIPE_EMAIL_TYPES}" \
    -J "01_Prepare_PAML_Sequences" \
    -o "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/01_Prepare_PAML_Sequences.stdout" \
    -e "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/01_Prepare_PAML_Sequences.stderr" \
    -N "${_PIPE_NNODES}" \
    -n "${_PIPE_NTASKS}" \
    -c "${_PIPE_CPUS_PER_TASK}" \
    -t "${_PIPE_WALLTIME}" \
    --mem-per-cpu "${_PIPE_MEM_PER_CPU}" \
    -p "${_PIPE_PARTITION}" \
    --export="_PIPE_SCRIPTS_FROM_GITHUB=${_PIPE_SCRIPTS_FROM_GITHUB},_PIPE_SCRATCH_DIR=${_PIPE_SCRATCH_DIR},_PIPE_RUN_NICKNAME=${_PIPE_RUN_NICKNAME},_PIPE_ALL_DATA=${_PIPE_ALL_DATA},_PIPE_FINAL_OUTPUT_DIR=${_PIPE_FINAL_OUTPUT_DIR},_PIPE_COHORT_MEMBERS=${_PIPE_COHORT_MEMBERS}" \
    "${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/Step_01_Prepare_PAML_Sequences.sh")
echo "Step 01: Prepare_PAML_Sequences has job ID ${STEP_01}" | tee -a "${_PIPE_EXEC_RECORD}"

# Step 02: Produce gene trees with RAxML
STEP_02=$(sbatch \
    --parsable \
    --kill-on-invalid-dep=yes \
    --dependency=afterok:"${STEP_01}" \
    --mail-type="${_PIPE_EMAIL_TYPES}" \
    -J "02_Make_Gene_Trees" \
    -o "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/02_Make_Gene_Trees.stdout" \
    -e "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/02_Make_Gene_Trees.stderr" \
    -N "${_PIPE_NNODES}" \
    -n "${_PIPE_NTASKS}" \
    -c "${_PIPE_CPUS_PER_TASK}" \
    -t "${_PIPE_WALLTIME}" \
    --mem-per-cpu "${_PIPE_MEM_PER_CPU}" \
    -p "${_PIPE_PARTITION}" \
    --export="_PIPE_FINAL_OUTPUT_DIR=${_PIPE_FINAL_OUTPUT_DIR}" \
    "${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/02_Make_Gene_Trees.sh")
echo "Step 02: Make_Gene_Trees has job ID ${STEP_02}" | tee -a "${_PIPE_EXEC_RECORD}"

# Step 03: Make PAML site model control files
STEP_03=$(sbatch \
    --parsable \
    --kill-on-invalid-dep=yes \
    --dependency=afterok:"${STEP_02}" \
    --mail-type="${_PIPE_EMAIL_TYPES}" \
    -J "03_Make_PAML_Control_Files" \
    -o "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/03_Make_PAML_Control_Files.stdout" \
    -e "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/03_Make_PAML_Control_Files.stderr" \
    -N "${_PIPE_NNODES}" \
    -n "${_PIPE_NTASKS}" \
    -c "${_PIPE_CPUS_PER_TASK}" \
    -t "${_PIPE_WALLTIME}" \
    --mem-per-cpu "${_PIPE_MEM_PER_CPU}" \
    -p "${_PIPE_PARTITION}" \
    --export="_PIPE_SCRIPTS_FROM_GITHUB=${_PIPE_SCRIPTS_FROM_GITHUB},_PIPE_FINAL_OUTPUT_DIR=${_PIPE_FINAL_OUTPUT_DIR}" \
    "${_PIPE_SCRIPTS_FROM_GITHUB}/Final_Pipeline_Scripts/03_Make_PAML_Control_Files.sh")
echo "Step 03: Make_PAML_Control_Files has job ID ${STEP_03}" | tee -a "${_PIPE_EXEC_RECORD}"

# #   Step 04: Broader metagenomics profiling
# STEP_04=$(sbatch \
#     --parsable \
#     --kill-on-invalid-dep=yes \
#     --dependency=afterok:"${STEP_03}" \
#     --mail-type="${_PIPE_EMAIL_TYPES}" \
#     -J "${_PIPE_NAME}.04_Metagenomics_Profiling" \
#     -o "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/04_Metagenomics_Profiling.stdout" \
#     -e "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/04_Metagenomics_Profiling.stderr" \
#     -N "${_PIPE_NNODES}" \
#     -n "${_PIPE_NTASKS}" \
#     -c "${_PIPE_CPUS_PER_TASK}" \
#     -t "${_PIPE_WALLTIME}" \
#     --mem-per-cpu "${_METAGENOME_MEM_PER_CPU}" \
#     -p "${_METAGENOME_PARTITION}" \
#     -A "${_PIPE_SLURM_ACCOUNT}" \
#     --export=_PIPE_WORKDIR="${_PIPE_WORKDIR}",_PIPE_INSTALL_DIR="${_PIPE_INSTALL_DIR}",_PIPE_CONDA_DIR="${_PIPE_CONDA_DIR}",_PIPE_RESULTS_DIR="${_PIPE_RESULTS_DIR}",_PIPE_NAME="${_PIPE_NAME}",_PIPE_METAGENOME_DB="${_PIPE_METAGENOME_DB}" \
#     "${_PIPE_INSTALL_DIR}/Pipeline_Scripts/04_Metagenomics_Profiling.sh")
# echo "Step 04: Metagenomics profiling has job ID ${STEP_04}" | tee -a "${_PIPE_EXEC_RECORD}"

# #   Step 05: Report generation. Note that we always use 1 CPU here.
# STEP_05=$(sbatch \
#     --parsable \
#     --kill-on-invalid-dep=yes \
#     --dependency=afterok:"${STEP_04}" \
#     --mail-type="${_PIPE_EMAIL_TYPES}" \
#     -J "${_PIPE_NAME}.05_Report" \
#     -o "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/05_Report.stdout" \
#     -e "${_PIPE_FINAL_OUTPUT_DIR}/Scheduler_Logs/05_Report.stderr" \
#     -N "${_PIPE_NNODES}" \
#     -n "${_PIPE_NTASKS}" \
#     -c "1" \
#     -t "${_PIPE_WALLTIME}" \
#     --mem-per-cpu "${_PIPE_MEM_PER_CPU}" \
#     -p "${_PIPE_PARTITION}" \
#     -A "${_PIPE_SLURM_ACCOUNT}" \
#     --export=_PIPE_WORKDIR="${_PIPE_WORKDIR}",_PIPE_INSTALL_DIR="${_PIPE_INSTALL_DIR}",_PIPE_CONDA_DIR="${_PIPE_CONDA_DIR}",_PIPE_RESULTS_DIR="${_PIPE_RESULTS_DIR}",_PIPE_NAME="${_PIPE_NAME}",_PIPE_FASTQ_DIR="${_PIPE_FASTQ_DIR}" \
#     "${_PIPE_INSTALL_DIR}/Pipeline_Scripts/05_Report.sh")
# echo "Step 05: Report generation has job ID ${STEP_05}" | tee -a "${_PIPE_EXEC_RECORD}"