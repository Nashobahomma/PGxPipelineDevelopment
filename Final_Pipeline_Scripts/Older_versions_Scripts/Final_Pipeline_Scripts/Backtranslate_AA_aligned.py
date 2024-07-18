#!/usr/bin/env python
"""
Authors: TK AN... Majority taken directly from code by T. Kono.
Date Modified: 2022-08-07

Back translate amino acid sequences to nucleotides for cyp-identified
orthogroups that have already been aligned. Requires Biopython and SAMtools.
Takes three arugments:
    1) Path to directory that contains CDS FASTA files
    2) Path to the aligned amino acid FASTA file
    3) Path to the Orthogroups.tsv file

Usage:

python /path/to/Backtranslate_AA_aligned.py /path/to/CDS_FASTA /path/to/OGXXXXXX_AA.fa /path/to/Orthogroups.tsv > /path/to/OGXXXXX_Backtranslated.fa 

"""

import sys
import subprocess
import os

#NOT DO THIS CONVENTION, SO IT IS MORE PORTABLE TO OTHER COMPUTING ENVIRONMENTS
#HOME_DIR = 
#OUTPUT_DIR = HOME_DIR +

try:
    cds_db_dir = sys.argv[1]
    og_to_backtrans = sys.argv[2]
    orthogroups_tsv = sys.argv[3]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    print('This script requires Biopython')
    exit(1)

def list_files(directory):
    """Get the FASTA files that are present in the supplied directory."""
    try:
        abpath = os.path.abspath(directory)
        filepaths = os.listdir(abpath)
        cds_dict = {}
        for fasta_file in filepaths:
            if fasta_file.endswith('.fasta'):
                species_name = fasta_file.split('.')[0]
                cds_dict[species_name] = os.path.join(abpath, fasta_file)
    except (OSError, IOError):
        print('The directory supplied is not readable, or does not exist!')
        exit(1)
    return cds_dict


def parse_og_tsv(og_table):
    """Parse the Orthogroups.tsv file and return dictionary of
    prot_id -> species. Use the column headers of the tsv file to identify the
    species."""
    # We will use the protein ID as the lookup name (called the "key" in Python
    # dictionaries). The referenced value (called the "value" in Python in
    # dictionaries) will be the species.
    key_dat = {}
    with open(og_table, 'rt') as f:
        for line_number, row in enumerate(f):
            if line_number == 0:
                header = row.strip().split('\t')
                # Remove the first column ("Orthogroup") because it is not of
                # interest for this step.
                species_names = header[1:]
            else:
                table_row = row.strip().split('\t')
                # When multiple proteins from the same species are in a single
                # orthogroup, Orthofinder reports them with a ', ' between
                # the protein IDs. To maintain the association between the
                # protein IDs and the species names, we will use the zip()
                # function.
                for pid_list, sp_name in zip(table_row[1:], species_names):
                    # Now, we will split the protein ID list
                    prot_ids = pid_list.split(', ')
                    for pid in prot_ids:
                        # push the protein ID -> species association into the
                        # dictionary
                        key_dat[pid] = sp_name
    return key_dat


#def fix_seqname(sname):
    #"""Don't need this cause it was needed for grass names quirkiness"""

def extract_cds(file_list, prot_id, protein_key):
    """Given the protein ID from the alignment and the file list, get the
    correct FASTA to extract from, then use samtools to fetch the CDS sequence.
    """
    # Get the species name associated with the protein ID
    species = protein_key.get(prot_id)
    # Get the filename associated with the species name
    cds_fasta = file_list.get(species)
    #   Then, buld the command line
    cmd = ['samtools', 'faidx', cds_fasta, prot_id]
    proc = subprocess.Popen(
        cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    out, err = proc.communicate()
    #   Then, process the output. We will remove the first line, since it is
    #   just the sequence ID. We will then put all the nucleotides into a single
    #   long string.
    lines = out.decode('utf-8').split('\n')
    cds = ''.join(lines[1:])
    return cds

def backtranslate(p_seq, n_seq):
    """Iterate through the aligned protein sequence, and replace the amino acids
    with codon triplets from the CDS file."""
    #   Keep track of the new sequence. Also keep track of which codon we are
    #   actually processing (gaps don't count)
    newseq = ''
    codon = 0
    for aa in p_seq:
        if aa == '-':
            #Altering from TK's code of inserting three dashes, we are inserting three ? as to be compatible with PAML input requirements
            newseq += '???'
        else:
            newseq += n_seq[codon*3:(codon*3) + 3]
            codon += 1
    return newseq


def generate_paml_name(orig_name, protein_species_key):
    """Generate a PAML-friendly sequence name. It will be in the format of 
        Gen.spe_PROTEIN.ID
    where 'Gen' is the first three letters of the genus name, 'spe' is the
    first three letters of the species epithet, and 'PROTEIN.ID' is the NCBI
    protein identifier."""
    # First, get the species name from the protein->species dictionary
    sp_name = protein_species_key.get(orig_name)
    # 2024-01-21: Change to a three letter genus and three letter species
    # epithet format instead. Also remove the "_counter" that NCBI appends
    # to the protein identifier.
    # We will convert this into the three-letter format
    genus, sp_epithet = sp_name.split('_')
    genus_species_code = genus[0:3] + '.' + sp_epithet[0:3]
    # Next, extract the protein ID from the original sequence name. Split the
    # sequence name on underscores, and take the three final parts.
    seq_parts = orig_name.split('_')
    # 2024-01-21: Also chop off the "_counter" bit that NCBI appends to the
    # protein ID.
    prot_id = '_'.join(seq_parts[-3:-1])
    # Join the three letter code to the protein ID and return the resulting
    # string.
    new_name = genus_species_code + '_' + prot_id
    return new_name


def main(db_dir, msa, og_tsv):
    """Main function."""
    #   Get the paths of the CDS files
    paths = list_files(db_dir)
    # Parse the Orthogroups.tsv file and store it as a dictionary
    prot_og_key = parse_og_tsv(og_tsv)
    #   Parse the alignment
    aln = SeqIO.parse(msa, 'fasta')
    for sequence in aln:
        #s, p = fix_seqname(sequence.id)
        cdsseq = extract_cds(paths, sequence.id, prot_og_key) 
        bt_seq = backtranslate(sequence.seq, cdsseq)
        # Added on 2022-08-07: new function to generate a PAML-friendly sequence name
        paml_name = generate_paml_name(sequence.id, prot_og_key)
        #   Print out the sequence in FASTA format
        print('>' + paml_name + '\n' + bt_seq)
    return


main(cds_db_dir, og_to_backtrans, orthogroups_tsv)