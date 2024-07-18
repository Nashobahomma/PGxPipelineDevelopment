#!/usr/bin/env python
"""
Authors: AN and TK
Date: 2022-07-13

Script that translates an unaligned nucleotide multi-FASTA file to amino acids.
This is used to translate orthogroup sequences into amino acids prior to
alignment with MAFFT. Alignment will happen with amino acid sequences so that
the codon structure of the genes is preserved for codon-level analyses that
will be done with PAML.

Requires Biopython.

Take one argument:
	1) Nucleotide FASTA file to translate.
"""

# First: import "standard library" modules
import sys

# Second: import non-standard and third-party modules.
#  Authors note: we will wrap the non-standard imports in a try/except block
#                so that we can print informative messages if someone does not
#                have the libraries available.
try:
	from Bio import SeqIO
	from Bio import SeqRecord
except ImportError:
	sys.stderr.write('This script requires the Biopython library.\n')
	sys.exit(1)

# Third: define the arguments that the script will process.
#  Authors note: we will wrap the argument definition into a try/except block
#                as well so that if a user does not give any arguments, we can
#                remind them how the script is supposed to be run.
try:
	nuc_fa = sys.argv[1]
except IndexError:
	sys.stderr.write(__doc__ + '\n')
	sys.exit(1)

# Now the boilerplate section is over. We can actually write the code that the
# script will need to do its work.

# Initialize an empty list to hold our translated sequences.
aa_seqs = []

# Open a 'read' handle to the nucleotide FASTA. The 'rt' argument in the
# open() function means "read" and "text" mode.
nuc_fa_handle = open(nuc_fa, 'rt')

# Then, actually read data out of the file handle. We will use the
# SeqIO.parse() function from Biopython. We will use this function
# because it is written specialized for reading FASTA files.
for nuc_seq in SeqIO.parse(nuc_fa_handle, 'fasta'):
	# Keep the sequence name and descriptions
	seq_name = nuc_seq.id
	seq_descr = nuc_seq.description
	# Translate the nucleotides to amino acids. We will use NCBI translation
	# table 1 (standard genetic code) to do the translation.
	aa_seq = nuc_seq.seq.translate(table=1)
	# Next, we will make a SeqRecord object that holds the translated sequence
	# and the original name. SeqRecord is a "class" of object defined as part of
	# the Biopython library.
	aa_seq_record = SeqRecord.SeqRecord(
		aa_seq,
		id=seq_name,
		description=seq_descr)
	# Append this SeqRecord object to the list defined on line 46
	aa_seqs.append(aa_seq_record)

# Now, we have translated all of the sequences in the input file to amino
# acids and are storing them in a list of SeqRecord objects. We will use the
# SeqIO.write() function in Biopython to export them as a properly-formatted
# FASTA file.
SeqIO.write(aa_seqs, sys.stdout, 'fasta')