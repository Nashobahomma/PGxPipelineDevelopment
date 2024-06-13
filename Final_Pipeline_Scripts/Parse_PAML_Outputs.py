#!/usr/bin/env python
"""Draft script to parse the PAML output for the CYP evolution project.
Date: 2023-09-24
Authors: AN and TK

Takes 2 arguments:
	1) Directory of PAML output files (.b)
	2) Orthogroup-CYP name CSV

Usage:
A one-line bash script that calls this script like so:
python /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/SevenSpeciesWorkflow/Step_07_ParseOutFiles/Draft_PAML_Output_Parser.py /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/SevenSpeciesWorkflow/Step_06_PAML_Runs /Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/SevenSpeciesWorkflow/TargetCYP_Sequences/CYPnames_Trans_Prot_with_OGs_7species.csv

Basically, this template:
python /path/to/Draft_PAML_Output_Parser.py /path/to/Step_06_PAML_Runs /path/to/CYPnames_Trans_Prot_with_OGs_7species.csv



"""

import sys
import os

try:
	paml_out_dir = sys.argv[1]
	cyp_name_csv = sys.argv[2]
except IndexError:
	sys.stderr.write(__doc__ + '\n')
	sys.exit(1)


def identify_paml_outputs(pdir):
	"""List the contents of the specified directory and return paths to the
	01278 and 8a output files in a dictionary. Isolate the orthogroup ID from
	the filename and use it as the key of the dictionary."""
	p_dict = {}
	for fname in os.listdir(pdir):
		# PAML output filenames end with '.out'
		if fname.endswith('.out'):
			# Figure out if we are handling a file for 01278 or 8a
			if fname.endswith('.01278.out'):
				mod = '01278'
			elif fname.endswith('.8a.out'):
				mod = '8a'
			# Isolate the orthogroup ID from the filename
			og_id = os.path.basename(fname).split('.')[0]
			# Build the full path to the output file and store it in the
			# dictionary
			full_path = os.path.join(
				os.path.abspath(os.path.expanduser(pdir)),
				fname)
			# Put the filename into the dictionary. Be sure to keep the 01278
			# and 8a files separate
			if og_id in p_dict:
				p_dict[og_id][mod] = full_path
			else:
				p_dict[og_id] = {mod: full_path}
	return p_dict


def parse_paml_output(paml_output: str) -> dict:
    """
    Parse the PAML output files. Return a dictionary with the model names
    as keys and the parameter estimates, lnL values, and significant BEB
    codons as the values.
    
    Args:
    paml_output (str): Path to the PAML output file to be parsed.
    
    Returns:
    dict: A dictionary with model names as keys and their parameter estimates,
          lnL values, and significant BEB codons as values.
    """
    # Initialize an empty dictionary to populate
    paml_dict = {}
    # Set the initial state for the current model and all parameters.
    # These variables are reset to the initial values each time a new model
    # section is encountered. This ensures each model is parsed independently.
    current_model = None
    beb_sites = None
    w, p, lnl, np = [], [], None, None

    def finalize_model() -> None:
        """
        Finalize the current model by computing dN/dS and storing the collected
        parameters in the dictionary.
        
        This step checks that the information parsed from the model contains
        values for lnl, w, and p. Models with missing values for these 
        parameters are not added to the final results dictionary. This means
        the directories from files ending in `01278` will not include Model 8a,
        and directories from files ending in `8a` will only include Model 8a.
        """
        if current_model and lnl and w and p:
            dn_ds = sum(float(weight) * float(site_class_omega) for weight, site_class_omega in zip(p, w))
            paml_dict[current_model] = {
                'w': w,
                'lnL': lnl,
                'np': np,
                'p': p,
                'dN/dS': str(dn_ds),
                'BEB.Significant': beb_sites
            }
    # Start parsing through the PAML file, line-by-line
    with open(paml_output, 'rt') as f:
        for line in f:
            stripped_line = line.strip()
            
            if stripped_line.startswith('ns = '):
                num_species = stripped_line.split()[2]
                paml_dict['Global'] = {'ns': num_species}
            if 'Model 0:' in stripped_line:
                finalize_model()
                current_model = 'Model 0'
                beb_sites = ['NA']
                w, p, lnl, np = [], [], None, None
            elif 'Model 1:' in stripped_line:
                finalize_model()
                current_model = 'Model 1'
                beb_sites = ['NA']
                w, p, lnl, np = [], [], None, None
            elif 'Model 2:' in stripped_line:
                finalize_model()
                current_model = 'Model 2'
                beb_sites = []
                w, p, lnl, np = [], [], None, None
            elif 'Model 7:' in stripped_line:
                finalize_model()
                current_model = 'Model 7'
                beb_sites = ['NA']
                w, p, lnl, np = [], [], None, None
            elif 'Model 8:' in stripped_line:
                finalize_model()
                current_model = 'Model 8'
                beb_sites = []
                w, p, lnl, np = [], [], None, None
            elif 'Model: One dN/dS' in stripped_line:
                finalize_model()
                current_model = 'Model 8a'
                beb_sites = ['NA']
                w, p, lnl, np = [], [], None, None
            if 'lnL' in stripped_line:
                lnl_pieces = stripped_line.split()
                np_piece = lnl_pieces[3]
                lnl = lnl_pieces[4]
                np = np_piece[:-2]
            if 'omega' in stripped_line and current_model == 'Model 0':
                w_pieces = stripped_line.split()
                w = [w_pieces[3]]
                p = ['1']
            if 'w:' in stripped_line and current_model != 'Model 0':
                w_pieces = stripped_line.split()
                w = w_pieces[1:]
            if 'p:' in stripped_line:
                p_pieces = stripped_line.split()
                p = p_pieces[1:]
            if 'Bayes Empirical Bayes' in stripped_line:
                for _ in range(5):
                    next(f)
                beb_text = next(f).strip()
                while not beb_text.startswith('The grid'):
                    beb_parts = beb_text.split()
                    if len(beb_parts) == 6 and '*' in beb_parts[2]:
                        beb_sites.append(beb_parts)
                    beb_text = next(f).strip()
            if 'Time used' in stripped_line:
                finalize_model()
                current_model = None  # Reset current_model after finishing parsing a model

    # Finalize the last model if not already done
    finalize_model()

    return paml_dict


def parse_cyp_names(c_names):
	"""Parse the OrhogroupID->CYP name CSV and return it as a dictionary. This
	is to include the CYP name in the collated PAML output table."""
	cyp_table = {}
	with open(c_names, 'rt') as f:
		for index, line in enumerate(f):
			if index == 0:
				continue
			csv_row = line.strip().split(',')
			# The columns are CYP name, Protein ID, Orthogroup ID]
			cyp_name = csv_row[0]
			prot_id = csv_row[1]
			og_id = csv_row[2]
			# Use the orthogroup ID as the key for the dictionary
			# 2024-03-23: Update the CYP name handling clode to handle multiple CYPs
			# in the same orthogroup. We will join the CYP names with underscores and
			# the protein IDs with semicolons (because the protein IDs have underscores)
			if og_id in cyp_table:
				# Extract the "old" CYP name and protein ID
				old_cyp_name = cyp_table[og_id][0]
				old_prot_id = cyp_table[og_id][1]
				# Append the "new" CYP name with an underscore and the "new"
				# protein ID with a semicolon.
				new_cyp_name = old_cyp_name + '_' + cyp_name
				new_prot_id = old_prot_id + ';' + prot_id
				# Push the new values into the dictionary
				cyp_table[og_id] = (new_cyp_name, new_prot_id)
			else:
				cyp_table[og_id] = (cyp_name, prot_id)
	return cyp_table


def main(paml_output_directory, cyp_names):
	"""Main function to handle the PAML parsing."""
	# List the contents of the PAML output directory and return paths to
	# the 01278 and 8a outputs for each orthogroup.
	paml_outputs_by_og = identify_paml_outputs(paml_output_directory)
	# Parse the CYP-OrthogroupID-ProteinID table
	cyp_name_dict = parse_cyp_names(cyp_names)
	# Print out a header for the table. It is long.
	header = [
		'Orthogroup.ID',
		'CYP.Name',
		'Human.Protein.ID',
		'Num.Species',
		'Model0.lnL',
		'Model0.dNdS',
		'Model0.np',
		'Model1.lnL',
		'Model1.dNdS',
		'Model1.np',
		'Model1.W',
		'Model1.P',
		'Model2.lnL',
		'Model2.dNdS',
		'Model2.np',
		'Model2.W',
		'Model2.P',
		'Model2.BEB.Position',
		'Model2.BEB.AminoAcid',
		'Model2.BEB.Probability',
		'Model2.BEB.W',
		'Model2.BEB.WSE',
		'Model7.lnL',
		'Model7.dNdS',
		'Model7.np',
		'Model7.W',
		'Model7.P',
		'Model8.lnL',
		'Model8.dNdS',
		'Model8.np',
		'Model8.W',
		'Model8.P',
		'Model8.BEB.Position',
		'Model8.BEB.AminoAcid',
		'Model8.BEB.Probability',
		'Model8.BEB.W',
		'Model8.BEB.WSE',
		'Model8a.lnL',
		'Model8a.dNdS',
		'Model8a.np',
		'Model8a.W',
		'Model8a.P'
		]
	# Woof.
	print(','.join(header))
	# Next, parse the PAML reports in a loop and print out the information
	# for each Orthogroup
	for og_id in paml_outputs_by_og:
		# Parse the PAML file that holds the results for models 01278 and
		# return them as a dictionary.
		sys.stderr.write(og_id + '\n')
		m_01278_parsed = parse_paml_output(paml_outputs_by_og[og_id]['01278'])
		m_8a_parsed = parse_paml_output(paml_outputs_by_og[og_id]['8a'])
		# Use the orthogroup ID to lookup the CYP name and protein ID
		cyp_name, protein_id = cyp_name_dict[og_id]
		# Look through the BEB.Significant key of the model 2 and model 8
		# outputs. If they are empty, then report 'NA' for the BEB columns.
		if m_01278_parsed['Model 2']['BEB.Significant'] == []:
			m_01278_parsed['Model 2']['BEB.Significant'] = [['NA']*6]
		if m_01278_parsed['Model 8']['BEB.Significant'] == []:
			m_01278_parsed['Model 8']['BEB.Significant'] = [['NA']*6]
		# Then print the body of the table
		to_print = [
			og_id,
			cyp_name,
			protein_id,
			m_01278_parsed['Global']['ns'],
			m_01278_parsed['Model 0']['lnL'],
			m_01278_parsed['Model 0']['dN/dS'],
			m_01278_parsed['Model 0']['np'],
			m_01278_parsed['Model 1']['lnL'],
			m_01278_parsed['Model 1']['dN/dS'],
			m_01278_parsed['Model 1']['np'],
			';'.join(m_01278_parsed['Model 1']['w']),
			';'.join(m_01278_parsed['Model 1']['p']),
			m_01278_parsed['Model 2']['lnL'],
			m_01278_parsed['Model 2']['dN/dS'],
			m_01278_parsed['Model 2']['np'],
			';'.join(m_01278_parsed['Model 2']['w']),
			';'.join(m_01278_parsed['Model 2']['p']),
			';'.join([i[0] for i in m_01278_parsed['Model 2']['BEB.Significant']]),
			';'.join([i[1] for i in m_01278_parsed['Model 2']['BEB.Significant']]),
			';'.join([i[2] for i in m_01278_parsed['Model 2']['BEB.Significant']]),
			';'.join([i[3] for i in m_01278_parsed['Model 2']['BEB.Significant']]),
			';'.join([i[5] for i in m_01278_parsed['Model 2']['BEB.Significant']]),
			m_01278_parsed['Model 7']['lnL'],
			m_01278_parsed['Model 7']['dN/dS'],
			m_01278_parsed['Model 7']['np'],
			';'.join(m_01278_parsed['Model 7']['w']),
			';'.join(m_01278_parsed['Model 7']['p']),
			m_01278_parsed['Model 8']['lnL'],
			m_01278_parsed['Model 8']['dN/dS'],
			m_01278_parsed['Model 8']['np'],
			';'.join(m_01278_parsed['Model 8']['w']),
			';'.join(m_01278_parsed['Model 8']['p']),
			';'.join([i[0] for i in m_01278_parsed['Model 8']['BEB.Significant']]),
			';'.join([i[1] for i in m_01278_parsed['Model 8']['BEB.Significant']]),
			';'.join([i[2] for i in m_01278_parsed['Model 8']['BEB.Significant']]),
			';'.join([i[3] for i in m_01278_parsed['Model 8']['BEB.Significant']]),
			';'.join([i[5] for i in m_01278_parsed['Model 8']['BEB.Significant']]),
			m_8a_parsed['Model 8a']['lnL'],
			m_8a_parsed['Model 8a']['dN/dS'],
			m_8a_parsed['Model 8a']['np'],
			';'.join(m_8a_parsed['Model 8a']['w']),
			';'.join(m_8a_parsed['Model 8a']['p'])
		]
		# WOOF
		print(','.join(to_print))
	return


main(paml_out_dir, cyp_name_csv)