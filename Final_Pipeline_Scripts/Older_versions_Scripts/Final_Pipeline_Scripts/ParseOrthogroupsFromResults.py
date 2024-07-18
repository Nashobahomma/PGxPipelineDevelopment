#!/usr/bin/env python
"""Identify and copy orthogroup FASTA files with human CYPs of interest into a
specified directory. Takes four arguments:
    1) CSV of CYPs to keep
    2) Orthogroups.tsv path
    3) Directory of Orthogroup sequences
    4) Destination directory

The CSV should have two columns:
    1: CYP name
    2: NCBI protein ID (for human) 

Usage: python /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/ParseOrthogroupsFromResults.py /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/CYPnames_Trans_Prot.csv /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/Dec18_Orthogroups.tsv > /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/CYP_Orthogroup_Table.csv                                    

inputs etc in readable english"""

# Edit: 2023-01-08
#       Added functionality to copy target orthogroup sequences into a destination
#       directory.

import sys
import pprint
import shutil
import os

try:
    cyp_table = sys.argv[1]
    og_table = sys.argv[2]
    og_results_dir = sys.argv[3]
    cyp_dest_dir = sys.argv[4]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def parse_cyp_table(ct):
    """Parse a CSV that holds the CYP-protein ID assocation. The first column
    should be the CYP name and the second column should be the NCBI protein ID.
    Returns a dictionary with the protein ID as the key name and the CYP name
    as the value."""
    try:
        fh = open(ct, 'rt')
        fh.close()
    except (FileNotFoundError, OSError):
        sys.stderr.write(
            'The CYP translation table does not exist, or is not readable.\n')
        sys.exit(1)
    proteins_of_interest = dict()
    with open(ct, 'rt') as f:
        for line_number, line in enumerate(f):
            if line_number == 0:
                continue
            cyp_name, prot_id = line.strip().split(',')
            proteins_of_interest[prot_id] = cyp_name
    return proteins_of_interest


def get_og_fasta_files(orthofinder_results_dir, og_list, target_cyp_dir):
    """Crawl through the orthofinder results directory to copy the target orthogroup
    FASTA files into a special destination directory, to be copied to GCS."""
    try:
        orthofinder_ab_path = os.path.abspath(os.path.expanduser(orthofinder_results_dir))
        filepaths = os.listdir(orthofinder_ab_path)
    except (OSError, IOError):
        print('The Orthofinder results directory supplied is not readable, or does not exist!')
        exit(1)
    try:
        target_ab_path = os.path.abspath(os.path.expanduser(target_cyp_dir))
    except (OSError, IOError):
        print('The target CYP directory supplied is not readable, or does not exist!')
        exit(1)
    for fasta_file in filepaths:
        if fasta_file.endswith('.fa'): #ORTHOGROUPS END IN .FA at present stage
            orthogroup_name = fasta_file.split('.')[0] #doesnt exactly matter, just taking the OG name not file extension or the period
            # If the OG name is in the list of OGs of interest (orthogroups that have CYP genes)
            if orthogroup_name in set(og_list.values()):
                # Copy the orthogroup sequence file into the destination dir
                # Use the shutil.copy2() function so that the file metadata, such as
                # acess and modifictaion timestamps, are preserved. File timestamps
                # are useful to distinguish when a given file was made and when a
                # particular analysis was performed.
                #   First: build the destination filename
                dest_fname = os.path.join(target_ab_path, orthogroup_name + '.fa')
                #   Second: build the source filename
                src_fname = os.path.join(orthofinder_ab_path, orthogroup_name + '.fa')
                #   Third: copy the file
                shutil.copy2(src_fname, dest_fname)
    return


def scan_ogs_for_cyps(c_dict, og_table):
    """Scan through the Orthogroups.tsv file to identify orthogroups that have
    human CYP genes of interest. Return a dictionary keyed on CYP name that
    has the orthogroup data as the value."""
    try:
        fh = open(og_table, 'rt')
        fh.close()
    except (FileNotFoundError, OSError):
        sys.stderr.write(
            'The Orthogroups.tsv file does not exist, or is not readable.\n')
        sys.exit(1)
    ogs_of_interest = dict()
    with open(og_table, 'rt') as f:
        for line_number, line in enumerate(f):
            if line_number == 0:
                # line is a string-type variable. string.strip() [strip is the
                # method] returns a copy of 'string' with whitespace removed
                # from the ends. string.split('\t') will "break" a string on
                # tab characters, and return a list.
                # Thus, line.strip().split('\t') returns a list of column
                # names from the TSV file.
                header = line.strip().split('\t')
                # Find out which column is called 'Homo_sapiens'
                try:
                    hs_column = header.index('Homo_sapiens')
                except ValueError:
                    sys.stderr.write(
                        'The Orthogroups.tsv file does not have a "Homo_sapens" column.')
                    sys.exit(1)
            else:
                # Explicitly strip only newlines from the Orthogroups.tsv file because
                # cases in which a species does not have proteins in a OG, it will
                # have just an empty column. string.strip() will remove these,
                # but they are meaningful, so we need to retain them.
                og_row = line.strip('\n').split('\t')
                # Orthogroup ID is the first column
                og_id = og_row[0]
                # Isolating the human protein IDs that are part of this orthogroup
                hum_prots = og_row[hs_column]
                # Split the comma+space separated string into a list
                hum_prots_list = hum_prots.split(', ')
                for hp in hum_prots_list:
                    # Isolate just the NCBI protein ID from the complicated name
                    # that the proteins have in the Orthogroups.tsv file:
                    # Orthogroups.tsv has entries like lcl|NC_000019.10_cds_NP_001268900.1_117358
                    # We need just the "NP_001268900.1" bit.
                    p_id = '_'.join(hp.split('_')[-3:-1])
                    # Now we can check if the protein ID is in the set of CYPs
                    # of interest
                    if p_id in c_dict:
                        # 2024-02-18: Edit to report all human CYPs in an orthogroup if
                        # there are multiple human CYPs in the orthogroup. This addresses
                        # the issue of 'NA' appearing in the CYP->OG CSV table.
                        cyp_of_interest = p_id
                        ogs_of_interest[cyp_of_interest] = og_id
                    else:
                        continue
    # Now we are done scanning the file. Return the dictionary
    return ogs_of_interest


def print_cyp_ogs(c_dict, co_dict):
    """Print a CSV that contains CYP name, NCBI protein ID, and orthogroup ID.
    We will also print a header."""
    csv_header = ['CYP_Name', 'Protein_ID', 'Orthogroup_ID']
    print(','.join(csv_header))
    # Iterate through the CYP-Protein dictionary (from the user-supplied CSV)
    for cyp_ncbi_id in c_dict:
        cyp_name = c_dict[cyp_ncbi_id]
        # When getting the orthogroup ID, we will use the dict.get() method
        # so that we can return a default value in the case that the CYP did
        # not appear in an orthogroup.
        cyp_og_id = co_dict.get(cyp_ncbi_id, 'NA')
        # Print the CSV row
        csv_row = [cyp_name, cyp_ncbi_id, cyp_og_id]
        print(','.join(csv_row))
    return


def main(c_table, o_table, o_dir, dest_dir):
    """Main function"""
    # Read the data from the CYP-Protein ID CSV into a dictionary
    cyp_dict = parse_cyp_table(c_table)
    # Scan through the orthogroups.tsv file and print the information for
    # those with CYPs of interest
    cyp_og_dict = scan_ogs_for_cyps(cyp_dict, o_table)
    # Copy the orthogroup sequences into the destination directory
    get_og_fasta_files(o_dir, cyp_og_dict, dest_dir)
    # Next, print a CSV (with header) of CYP name, protein ID, and orthogroup ID
    print_cyp_ogs(cyp_dict, cyp_og_dict)
    return



main(cyp_table, og_table, og_results_dir, cyp_dest_dir)