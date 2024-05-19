#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem-per-cpu=4gb
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --time=16:00:00 
#SBATCH --mail-user=amber.nashoba@cuanschutz.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Make_Gene_Trees_21Apr2024_%j.err
#SBATCH -o Make_Gene_Trees_21Apr2024_%j.out

# This script generates gene trees for all target ortholog/cyp multiple sequence alignments
# (msa) files. These gene trees are one of the three required elements to run PAML analysis.
# On the PAML help page, PAML author Yang recommends using the best tree for the gene of # interest. 
# Requires install of RAXML into $PATH (i.e. needs to be installed into a location in which it can
# be called directly from a terminal window).
# This code is meant to be run from you local work station using previously generated PAML input
# sequence files.
# usage: bash /PATH/TO/Make_Gene_Trees.sh

# JOB COMMANDS
module load anaconda
conda activate CYPevol
# create a symbolic link to the data in the scratch directory and make sure any place files are read read
# from the scratch path. Or copy the data to the scratch directory and read from scratch.
cp -r /projects/anashoba@xsede.org/ForAlpine/Step_03_PAML_Seq_Inputs /scratch/alpine/$USER
cp -r /projects/anashoba@xsede.org/ForAlpine/Step_04_PAML_Gene_Trees /scratch/alpine/$USER


#RAXML_INPUT="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_03_PAML_Seq_Inputs"
#TREE_OUTPUT="/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_04_PAML_Gene_Trees"
RAXML_INPUT="/scratch/alpine/anashoba@xsede.org/Step_03_PAML_Seq_Inputs"
TREE_OUTPUT="/projects/anashoba@xsede.org/ForAlpine/Step_04_PAML_Gene_Trees"

mkdir -p "${TREE_OUTPUT}" #-p means if the folder exists, move on and don't throw an error.

for gene_alignment in $(find "${RAXML_INPUT}" -mindepth 1 -maxdepth 1 -type f -name '*.fa')
	do 
		orthogroup_id="$(basename ${gene_alignment} | cut -d '_' -f 1)"
		raxml-ng --all --msa "${gene_alignment}" --msa-format FASTA --data-type DNA --threads "${SLURM_CPUS_PER_TASK}" \
	    --model GTR+G --bs-trees 10000 --prefix "${TREE_OUTPUT}/${orthogroup_id}"
	done


