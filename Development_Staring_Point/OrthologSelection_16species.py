#!/usr/bin/env python
"""
Select a "representative" orthlogue from each species of an aligned (and
backtranslated) orthogroup. The method that will be used for choosing a
representative is the sequence with the fewest number of gaps will be
retained. (Recall that because the sequences are aligned, they will all be
the same length.) In the case of a tie, the first sequence in lexicographic
sort order will be chosen. Requires Biopython, and take one argument:
	1) Backtranslated orthogroup FASTA for an individual gene
	2) CYP-ProteinID-Orthogroup CSV

Usage:
python OrthologSelection.py /path/to/OGXXXXX_Backtranslated.fa path/to/ProteinID_JulyOG.csv > /path/to/OGXXXXX_RepOrthologues.fa


2024-02-18: Update to the selection procedure to perform the following:
	1: Calculate ungapped length for each non-human gene sequence
	2: Retain the longest sequence from each non-human primate (as calcualted by ungapped length)
	3: Retain all human gene sequences
	4: Apply a column-wise filter such that columns with more than 50% gaps (? character)
		are removed.

In applying this procedure, we are specifically setting up for analysis of *human* sequences.
Note that we are choosing to discard gene sequences from non-human primates because we are
primarily interested in analyzing the genes as they are characterized in humans.

The column-wise gap filter is applied because "gappy" columns take a long time to
process in PAML and do not provide much information.

"""

import sys
import subprocess
import os
import pprint

try:
   from Bio import SeqIO
except ImportError:
   print('This script requires Biopython')
   sys.exit(1)

try:
	og_fa_in = sys.argv[1]
	cyp_protid_og_csv = sys.argv[2]
except IndexError:
	sys.stderr.write(__doc__+ '\n')
	sys.exit(1)

def generate_fasta_dictionary(og_fasta):
	"""Read a FASTA file containing a backtranslated and aligned orthogroup sequence. It will return a
	dictionary of the form
		
			'sp1_id': [sp1_gene1, sp1_gene2, ...],
			'sp2_id': [sp2_gene1, sp2_gene2, ...]
		}
	where 'spX_id' is the 6-character species identifier for species X, and 'spX_geneY' is the 'SeqRecord'
	object for the actual FASTA record of a single gene in species X.
	EDIT: 12 Feb 2024 - with the addition of species. now using a gen.spe six letter plus . code. try for 7 characters 0:6 ..."""
	og_dict = dict()
	for og_seq in SeqIO.parse(og_fasta, 'fasta'):
		seq_name = og_seq.id 
		species_name = seq_name[0:7]
		if species_name in og_dict: 
			og_dict[species_name].append(og_seq)
		else:  
			og_dict[species_name] = [og_seq]
	return og_dict


def get_orthogroupID(pathTo_fasta):
	"""path to a fasta file and return the orthogroup id that the input path corresponds to"""
	#basename of the input file:
	#>>> import os
	#>>> pth = '/home/tomkono/whatever/data.txt'
	#>>> os.path.basename(pth) 
	fasta_basename = os.path.basename(pathTo_fasta) #name of orthogroup without directory path 
	og_ID = fasta_basename.split('_')[0]
	return og_ID

def generate_gene_lookup(csvPath):
	"""return a dictionary of orthogroup-human protein association pairs"""
	lookup = dict()
	with open(csvPath, 'rt') as file_handle: #r=read, t=text
		for line_number, csv_row in enumerate(file_handle): #for each line in the file...
			if line_number == 0:
				continue #skip the first line because it contains only column names
			cyp_name, prot_id, orthogroup_id = csv_row.strip().split(',') #.strip gets rid of the 'newline' at end of csv row
			# Edit 2024-02-18: key on the protein IDs in the CSV, rather than the orthogroup IDs.
			# we know that protein IDs will be unique, and that there can be multiple proteins of
			# interest in a single orthogroup.
			lookup[prot_id] = orthogroup_id #using orthogroup id as the key because will need to lookup by orthogroup id
	return lookup

def pick_representative_ortholog(full_og_dict, og_ID, og_human_lookup):
	"""This function will parse the orthologs in each species group and 1)calculate the length of the
	sequence 2)calculate how many gap characters (????) are in the sequence 3)return the ortholog per
	species with the fewest gap characters."""
	# Uncomment this to see the format of the NCBI Protein ID->Orthogroup dictionary
	# pprint.pprint(og_human_lookup)
	representative_og_dict = dict()
	for species_group in full_og_dict:
		if species_group == 'Hom.sap': #look up reference gene-protein ID EDIT 12 feb 2024 to longer species name convention
			# Edit 2024-02-18: Keep all human genes that were identified in the CSV
			representative_og_dict['Hom.sap'] = []
			for human_gene in full_og_dict['Hom.sap']:
				# Extract the NCBI protein ID from the FASTA sequence identifier. The first
				# eight characters are the genus, full stop, species, and underscore.
				human_protid = human_gene.id[8:]
				if human_protid in og_human_lookup:
					representative_og_dict['Hom.sap'].append(human_gene)
		else:
			species_genelist = full_og_dict[species_group]
			seq_name_gap_count = []
			for sp_gene in species_genelist:
				seq_name = sp_gene.id
				gap_count = sp_gene.seq.count('?')
				seq_name_gap_count.append((seq_name, gap_count))
			seq_name_gap_count.sort(key=lambda x: x[1], reverse=False) #sort from low gap to high gap count
			representative_seq_id = seq_name_gap_count[0][0] #pick first pair and first member of the pair
			for sp_gene in species_genelist:
				if sp_gene.id == representative_seq_id:
					representative_og_dict[sp_gene.id] = sp_gene
	# Uncomment this to see the format of the representative orthogroup dictionary
	# pprint.pprint(representative_og_dict)
	return representative_og_dict


def column_filter(repr_dict, max_prop_gap=0.5):
	"""Apply a column-wise gapping filter to a dictionary of aligned nucleotide
	sequences. Columns with more than 'max_prop_gap' gap characters will be
	trimmed."""
	# Start to unpack the dictionary of represenatative sequences. The nonhuman
	# sequences are SeqRecord objects and the human sequences are stored as a list.
	sequence_matrix = dict()
	for sp in repr_dict:
		if sp == 'Hom.sap':
			for prot in repr_dict[sp]:
				sequence_matrix[prot.id] = str(prot.seq)
		else:
			sequence_matrix[sp] = str(repr_dict[sp].seq)
	# Now we have a dictionary of strings, where the key is the sequence name
	# (which includes genus, species, and NCBI protein ID) and the value is the
	# nucleotide sequence. We now have to apply the column-wise filter.
	# Let's start the filtering procedure by "transposing" the sequence matrix so
	# we can iterate over columns, rather than rows.
	t_seq_matrix = zip(*sequence_matrix.values())
	# Then, let's make a list of alignment columns that pass our gapping filter. We
	# will do this by iterating over the transposed sequence matrix and keeping a record
	# of column numbers that pass.
	keep_columns = []
	for col_number, aln_column in enumerate(t_seq_matrix):
		# Calculate the number of gaps (? characters)
		gap_count = aln_column.count('?')
		gap_prop = gap_count / len(aln_column)
		if gap_prop <= max_prop_gap:
			keep_columns.append(col_number)
	# Then, start a new dictionary for the filtered alignment. Keep only the columns
	# that passed our filter in the previous block of statements. We will refer to the
	# un-transposed alignment matrix now.
	filtered_seq_matrix = dict()
	for sequence in sequence_matrix:
		# Filter the aligned nucleotide row based on the columns that passed the gap filter.
		flt_row = ''
		for col_number, nucleotide in enumerate(sequence_matrix[sequence]):
			if col_number in keep_columns:
				flt_row += nucleotide
		# Push the filtered sequence into the new dictionary
		filtered_seq_matrix[sequence] = flt_row
	# Return the filtered alignment dictionary
	return filtered_seq_matrix


def mask_stop_codons(og_seq_dictionary):
	"""Screen for in-frame STOP codons and replace them with ambiguity codes
	('???'). We need this function because NCBI sometimes serves CDS reference
	sequences with internal STOP codons."""
	# Define the set of mammalian STOP codons
	stop_codons = set(['TAA', 'TAG', 'TGA'])
	masked_sequences = dict()
	# Iterate through each sequence
	for sequence in og_seq_dictionary:
		# Start an entry in the masked sequence dictionary
		masked_sequences[sequence] = ''
		# Then, iterate through the nucleotides in each sequence
		#	Annotation for the following statements:
		#		Iterate over the nucleotides with a counter starting at 1 (instead of 0) 
		#		If the counter, when divided by 3 has a remainder of 0 (i.e., we are on a multiple of 3)
		#				N.B.: This is called the "modulus" operator. X % Y = "X modulus Y" = "remainder when X is divided by Y"
		#			Then, take the current nucleotide and the two before the current and treat them as a codon
		#			If the codon is in the set of STOP codons, replace them with '???'
		#			Else, leave them as-is
		#			Put the codon onto the end of a growing string that holds the masked sequence
		#		Return the masked sequences
		for nuc_index, base in enumerate(og_seq_dictionary[sequence], start=1):
			if nuc_index % 3 == 0:
				codon = og_seq_dictionary[sequence][nuc_index-3:nuc_index]
				if codon in stop_codons:
					masked_sequences[sequence] += '???'
				else:
					masked_sequences[sequence] += codon
	return masked_sequences




def print_best_orthologs(dict_of_representative_orthologs):
	"""print a fasta file to standard output"""
	for sequence in dict_of_representative_orthologs:
		print('>' + sequence)
		print(dict_of_representative_orthologs[sequence])

def main(og_file, cyp_protein_orthogroup_table):
	"""function doc string: highlevel function that will carry out the steps of the algorithm described 
	in the module doc string above. The main function is pretty short and calls other functions defined
	in the script. What is happening in this script can be seen in main"""
	#read the orthogroup alignment into memory. Stored as as dictionary.
	og_ID = get_orthogroupID(og_file) # function to extract orthogroup id from the input fasta file name:
	og_human_lookup = generate_gene_lookup(cyp_protein_orthogroup_table)
	og_species_dict = generate_fasta_dictionary(og_file) #pseudocode step 1
	og_rep_per_species = pick_representative_ortholog(og_species_dict, og_ID, og_human_lookup) # pseudocode actions for step 2
	og_col_filtered = column_filter(og_rep_per_species)
	# 2024-03-17: Make a new function that screens for in-frame STOP codons and replaces them with '???'
	og_col_stop_filtered = mask_stop_codons(og_col_filtered)
	print_best_orthologs(og_col_stop_filtered)


main(og_fa_in, cyp_protid_og_csv)
