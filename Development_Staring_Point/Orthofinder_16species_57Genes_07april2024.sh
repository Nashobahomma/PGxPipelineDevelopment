#!/bin/bash

# JOB RESOURCES
###############

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem-per-cpu=4gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=16:00:00 
#SBATCH --mail-user=amber.nashoba@cuanschutz.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Orthofinder_07Apr2024_%j.err
#SBATCH -o Orthofinder_07Apr2024_%j.out

# Usage:
# sbatch /path/to/Orthofinder_16species_57Genes_07april2024.sh 
#	e.g., sbatch /projects/anashoba@xsede.org/ForAlpine/Scripts_March2024_4Alpine/Scripts_16species_57Genes_versions_January2024_OFruns/Orthofinder_16species_57Genes_07april2024.sh
# Note that Slurm will write the '.err' and '.out' files into the working directory of
# the shell.


# JOB COMMANDS
module load anaconda
conda activate CYPevol
# create a symbolic link to the data in the scratch directory and make sure any place files are read read
# from the scratch path. Or copy the data to the scratch directory and read from scratch.
cp -r /projects/anashoba@xsede.org/ForAlpine/CDS_March2024_4Alpine /scratch/alpine/$USER
cp /projects/anashoba@xsede.org/ForAlpine/CYP57genesHumanProteinIDs.csv /scratch/alpine/$USER

# Use Bash "strict mode"
#       -e: Exit on first error, rather than try to run the whole script
#       -u: Make undefined variables errors, rather than using empty strings as default value
#       -o pipefail: when one command in a pipe fails, treat all commands as having failed
set -euo pipefail

# Update 2023-03-26:
# Define the current run date as a variable that we can reference later in the script, so that if a run
# crosses the midnight time, then the dates will be consistent. This is for making sure that all ouput
# files that are named with a date stamp have consistent dates.
TODAY=$(date +'%F')
# Update 2024-02-18:
# Add a "tag" or a "nickname" variable to help identify the output files in Google Cloud
# bucket. This is to help distinguish orthofinder runs with different input data sets.
# i.e., output files will have this "NICKNAME" variable text in their names, like
# 2024-02-18_Orthogroups_16Species_57Genes.tsv
NICKNAME="16Species_57Genes"

# 2023-02-12: Update to separate outputs into three distinct parts:
#				1: Updated CSV with CYP name, NCBI protein ID, Orthogroup ID
#				2: Tarball with updated CSV and FASTA files of orthogroups with CYPs of interest
#				3: Tarball with archived orthofinder results
# 2022-12-18: Update to prepending a RAxML directory to the $PATH variable, as
# a dependency for Orthofinder.
#export PATH=/home/shared/amber.nashoba_PGx_Pipeline_Software/standard-RAxML-master:${PATH}

# 12/2022: This run is a trial run with time= geared for time to run only 3 species.

# Define paths/URLs to input data. The GCS names written into this example are for
# Amber Nashoba's test runs. Change the "Assemblies/" and "Orthofinder_Results/" parts of
# the GCS_INPUT_LOC and GCS_OUTPUT_LOC variables to match your folder names
#       URL to the Google Cloud Server bucket that has your CDS fasta files as 
#   directed in the Pipeline Workflow Document
GCS_INPUT_LOC="/scratch/alpine/anashoba@xsede.org/CDS_March2024_4Alpine/"

#	2024-01-21: This folder now has 16 FASTA files, so there are 16 species as
#	input for this pipeline.
#       Name of the folder in the bucket that has the CDS fast files. We have to
#       define this as a separate variable because it will be reused as the name
#       of a directory on the HPC cluster.
CDS_FOLDER="/scratch/alpine/anashoba@xsede.org/CDS_March2024_4Alpine/CDS_March2024_16sp"

#
#		URL to Google Cloud Server bucket for targeted orthogroup files (e.g., CYPs of interest)
mkdir -p /scratch/alpine/anashoba@xsede.org/Targeted_OG_Sequences
TARGET_OUTPUT_LOC="/scratch/alpine/anashoba@xsede.org/Targeted_OG_Sequences"

#       Edit on 2023-02-12: Define a name for the directory into which we will copy the targeted
#		OGs.
TARGET_CYP_DIR_BASENAME="TargetCYP_Sequences"

# Define paths to programs and resources on the HPC cluster
#       Path to the Orthofinder Singularity image.
# ORTHOFINDER="/home/shared/amber.nashoba_PGx_Pipeline_Software/orthofinder_tw_20220623.sif"
#       Edit on 2023-01-08: Path to ParseOrthogroupsFromResults.py, to isolate orthogroups
#       of interest based on having a CYP.
PARSE_ORTHOGROUPS_PY="/projects/anashoba@xsede.org/ForAlpine/Scripts_March2024_4Alpine/Scripts_16species_57Genes_versions_January2024_OFruns/ParseOrthogroupsFromResults.py"
#       Edit on 2023-01-08: Path to the two-column CSV (CYP name and human protein ID) for
#       identifying orthogroups of interest. Each user will have a personal CSV with their
#       genes of interest. This file is NOT in the group's shared directory.
CYP_PROT_CSV="/scratch/alpine/anashoba@xsede.org/CYP57genesHumanProteinIDs.csv"
# Fetch input data from GCS into a temporary directory on the HPC cluster
#       Write a progress message to the .out file
echo "$(date +'%F %T'): Making temp directory for analysis"
#       Use the mktemp comand to return a directory path in the system temp space
cd /scratch/alpine/anashoba@xsede.org
d=$(mktemp -d)
cd "${d}"


# Run Orthofinder
#	Define a output folder name for Orthofinder based on the date that the run
#	is performed.
echo "$(date +'%F %T'): Making output directory for Orthofinder"
HPC_ORTHOFINDER_OUT="$(date +'%F')_Orthofinder_Results_B2_${NICKNAME}"
#	Remove the output folder if already exists, to prevent "cross-contamination" of
#	Orthofinder analyses.
rm -rf "./${HPC_ORTHOFINDER_OUT}"
echo "$(date +'%F %T'): Starting Orthofinder"
orthofinder -f "${CDS_FOLDER}" -t "${SLURM_CPUS_PER_TASK}" -a "${SLURM_CPUS_PER_TASK}" -d -M "msa" -os -S "blast_nucl" -A "mafft" -T "raxml" -z -o "${HPC_ORTHOFINDER_OUT}" 
echo "$(date +'%F %T'): Finished Orthofinder"

# Edit on 2023-01-08: Add a step to run ParseOrthogroupsFromResults.py to isolate the
# orthogroups of interest, based on having a CYP.
echo "$(date +'%F %T'): Running ParseOrthogroupsFromResults.py to isolate OGs of interest"
# Use a find command to get the full path to the Orthogroups.tsv file in the Orthofinder results directory
# This is because Orthofinder names its output based on the date it is run, which will necessarily change each run.
OG_TSV=$(find "${HPC_ORTHOFINDER_OUT}" -type f -name 'Orthogroups.tsv')
echo "$(date +'%F %T'): Using ${OG_TSV} as the path for Orthogroups.tsv"
# Use a find command to get the path to the Orthogroup_Sequences directory in the Orthofinder results directory.
OG_SEQS_DIR=$(find "${HPC_ORTHOFINDER_OUT}" -type d -name 'Orthogroup_Sequences')
echo "$(date +'%F %T'): Using ${OG_SEQS_DIR} as the path for Orthogroup_Sequences"
# Make the full path to the destination directory for target CYPs. Use a find command to get the name of the
# date-specific results directory that is made by Orthofinder.
ORTHOFINDER_DATE_DIR=$(find "${HPC_ORTHOFINDER_OUT}" -mindepth 1 -maxdepth 1 -type d -name 'Results_*')
echo "$(date +'%F %T'): Using ${ORTHOFINDER_DATE_DIR} as the base path for Orthofinder results"
TARGET_CYP_DIR_FULLPATH="${d}/${TARGET_CYP_DIR_BASENAME}"
echo "$(date +'%F %T'): Using ${TARGET_CYP_DIR_FULLPATH} as the path for saving target CYP-containing orthogroups"
# Make the target directory
mkdir -p "${TARGET_CYP_DIR_FULLPATH}"
echo "$(date +'%F %T'): running 'python ${PARSE_ORTHOGROUPS_PY} ${CYP_PROT_CSV} ${OG_TSV} ${OG_SEQS_DIR} ${TARGET_CYP_DIR_FULLPATH} > ${TARGET_CYP_DIR_FULLPATH}/CYPnames_Trans_Prot_with_OGs_${NICKNAME}.csv'"
python "${PARSE_ORTHOGROUPS_PY}" "${CYP_PROT_CSV}" "${OG_TSV}" "${OG_SEQS_DIR}" "${TARGET_CYP_DIR_FULLPATH}" > "${TARGET_CYP_DIR_FULLPATH}/CYPnames_Trans_Prot_with_OGs_${NICKNAME}.csv"
echo "$(date +'%F %T'): Finished running ParseOrthogroupsFromResults.py to isolate OGs of interest"

# 2024-02-25
# Update: we are going to put the Orthogroups TSV into the same directory as the targeted
# CYP dir and CSV. This is to simplify the output transfer logistics.
	# Update on 2023-02-36
	# Rename the Orthogroups.tsv output file to include the run date and pushing to GCS.
# echo "$(date +'%F %T'): Renmaing the Orthogroups.tsv file to include the date and copying to GCS"
echo "$(date +'%F %T'): Copying the Orthogroups.tsv file into the targeted CYP directory"
# Update on 2024-02-25: Add the "nickname" (defined on line 27) to the TSV output to
# avoid clobbering outputs when multiple Orthofinder scripts are run on the same date.
cp "${OG_TSV}" "${TARGET_CYP_DIR_FULLPATH}/${TODAY}_Orthogroups_AllOGs_${NICKNAME}.tsv"
# 2024-02-25: Uncomment this command to send the Orthogroups TSV to the Google bucket
# gsutil -m cp -r "${TODAY}_Orthogroups_AllOGs_${NICKNAME}.tsv" "${TARGET_OUTPUT_LOC}"
# echo "Done pushing Orthogroups.tsv to GCS"

# Udpate on 2023-02-12
# Tar-gzip the targeted OG directory
echo "$(date +'%F %T'): Tar-gzipping the targeted OG directory"
TARGET_TARBALL="${TODAY}_Targeted_OGs_${NICKNAME}.tar.gz"
tar -cf - "./${TARGET_CYP_DIR_BASENAME}" | pigz -p "${SLURM_CPUS_PER_TASK}" -c > "${TARGET_TARBALL}"
echo "$(date +'%F %T'): Done tar-gzipping the targeted OG directory"

# Push the tarball of the target OGs to GCS
echo "$(date +'%F %T'): Pushing target OG tarball to targeted OG folder in GCS"
cp -r "${TARGET_TARBALL}" "/projects/anashoba@xsede.org"
echo "$(date +'%F %T'): Done pushing target OG tarball to targeted OG folder in GCS"

# Remove the temp directory that we used for analysis. Not strictly necessary,
# but nice for multiuser systems.
cd
rm -rf "${d}"
echo "$(date +'%F %T'): Done"