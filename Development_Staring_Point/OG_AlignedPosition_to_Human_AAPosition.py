#!/usr/bin/env python
"""Simple script to calculate a human codon position from a query aligned
amino acid position (e.g., one that is reported by PAML). Requires Biopython.
Takes two arguments:
	1) Alignment FASTA that was used for PAML
	2) Aligned codon position to convert

Example usage:
python /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/SeqFile_GenerationScripts_4paml/OG_AlignedPosition_to_Human_AAPosition.py /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/SevenSpeciesWorkflow/Step_03_PAML_Seq_Inputs/OG0016452_RepOrthologues.fa 54"""

import sys

try:
	aln_in = sys.argv[1]
	aligned_pos = int(sys.argv[2])
except IndexError:
	sys.stderr.write(__doc__ + '\n')
	sys.exit(1)
except ValueError:
	sys.stderr.write('The aligned codon position should be an integer.\n')
	sys.exit(1)

try:
	from Bio import SeqIO
	from Bio import Seq
except ImportError:
	sys.stderr.write('This script requires Biopython.\n')
	sys.exit(1)


def main(og_aln, paml_position):
	"""Main function."""
	# First, open a handle to the FASTA for reading
	fa_handle = open(og_aln, 'rt')
	# Then, open a FASTA parser that lets us iterate through the aligned
	# sequences, one at a time
	fa_parser = SeqIO.parse(fa_handle, 'fasta')
	# Iterate through the parser, and identify which sequence is the human
	# sequence. The way we built the pipeline, the human sequence will always
	# start with 'H.sa'
	for aln_seq in fa_parser:
		if aln_seq.id.startswith('H.sa'):
			# Start a counter for the human sequence. This counter will only
			# increment if we see a non-gap character in the human sequence.
			non_gap_counter = 0
			# Also keep track of the codon sequence; we will use to translate
			# to amino acid so we can report the human reference amino acid
			# state.
			codon_seq = ''
			# And keep track of the human nucleotide positions as well.
			hsa_nuc_pos = []
			# Use enumerate() here to keep track of where we are in the aligned
			# sequence (will keep track of gaps to match up with PAML position).
			for aln_column, nucleotide in enumerate(aln_seq.seq):
				# Skip gaps (? or -)
				if nucleotide in ('?', '-'):
					continue
				# Then, convert the non_gap_counter (which is keeping track of
				# nucleotides) to amino acid coordinates. Divide by 3, and
				# discard the remainder. Add 1 because the first codon is at
				# position 1, not 0.
				hsa_aa_pos = int(non_gap_counter / 3) + 1
				# Calculate the aligned (non human) codon position.
				aln_aa_pos = int(aln_column / 3) + 1
				# Then, if our hsa_aa_pos is the query position, then add the
				# current nucleotide to the codon we are building
				if aln_aa_pos == paml_position:
					codon_seq = codon_seq + nucleotide
					hsa_variant_pos = hsa_aa_pos
					hsa_nuc_pos.append(str(non_gap_counter))
				# Increment the counter for the nucleotide that we just
				# processed.
				non_gap_counter += 1
			# Print out the codon, its translation, and the position		
			codon_trans = Seq.Seq(codon_seq).translate()
			print('Human amino acid position:', hsa_variant_pos)
			print('Human reference amino acid:', codon_trans)
			print('Human reference codon:', codon_seq)
			print('Human reference CDS position (nucleotides):', ', '.join(hsa_nuc_pos))
	return


main(aln_in, aligned_pos)