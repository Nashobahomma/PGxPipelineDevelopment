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


def parse_paml_output(paml_output):
	"""Parse the PAML output files. Return a dictionary with the model names
	as keys and the parameter estimates, lnL values, and significant BEB
	codons as the values."""
	# Start an empty dictionary to populate
	paml_dict = {}
	# Start parsing through the PAML file, line-by-line
	with open(paml_output, 'rt') as f:
		for line in f:
			# We will extract the "ns" value: the number of species in the alignment
			if line.startswith('ns = '):
				# We will rely on the default behavior of Python's string.split() method,
				# which is to split on any whitespace. We will take the third element of
				# the resulting list.
				num_species = line.strip().split()[2]
				# Put this number into the dictionary
				paml_dict['Global'] = {'ns': num_species}
				# If we want to inclue the length of the alignment, then we will do it
				# here, too.
			# We will take advantage of the fact that the PAML output file is
			# organized into coherent sections. That is, we will not find
			# outputs for Model 1 in the section for Model 0. BEB is only
			# meaningful for models 2 and 8 (those with positively selected
			# sites in the model).
			if line.startswith('Model 0:'):
				current_model = 'Model 0'
				beb_sites = ['NA']
			elif line.startswith('Model 1:'):
				current_model = 'Model 1'
				beb_sites = ['NA']
			elif line.startswith('Model 2:'):
				current_model = 'Model 2'
				beb_sites = []
			elif line.startswith('Model 7:'):
				current_model = 'Model 7'
				beb_sites = ['NA']
			elif line.startswith('Model 8:'):
				current_model = 'Model 8'
				beb_sites = []
			elif line.startswith('Model: One dN/dS'):
				# This line is to identify model 8a, which was run with a separate
				# control file, so it has a slightly different format from the others.
				current_model = 'Model 8a'
				beb_sites = ['NA]']
			# Then, the lnL values and np (we think this is the number of parameters)
			# values, and the p/w values will be extracted in the same way for each
			# model.
			if line.startswith('lnL'):
				# Split the string on whitespace. The fourth item will have the
				# np value, and the fifth will have the lnL value.
				lnl_pieces = line.strip().split()
				np_piece = lnl_pieces[3]
				lnl = lnl_pieces[4]
				# Note that the np value will have a trailing '):' after the
				# number. So, we will exclude the final two characters.
				np = np_piece[:-2]
			# Model 0 does not have 'p' nor 'w' output lines. But, it does have
			# a line that starts with 'omega'
			if line.startswith('omega'):
				# We will name this parameter 'w' in our script to make it
				# consistent with the outputs from the other models
				w_pieces = line.strip().split()
				# Because model 0 reports one ratio, we must store it in a list of
				# length 1 to make it a consistent object type with the outputs from
				# the models that report multiple ratios
				w = [w_pieces[3]]
				# This is a "dummy" value for p (proportion of the gene with a given
				# dN/dS estimate) for model 0
				p = ['1']
			# The 1, 2, 7, and 8 models use 'w' and 'p' to report dN/dS and
			# the proportion of the gene with specific dN/dS ratios
			if line.startswith('w:'):
				w_pieces = line.strip().split()
				w = w_pieces[1:]
			if line.startswith('p:'):
				p_pieces = line.strip().split()
				p = p_pieces[1:]
			# If the line starts with 'Bayes Empirical Bayes' then we are in the BEB
			# section. This is only present in models 2 and 8 (positive selection).
			if line.startswith('Bayes Empirical Bayes'):
				# The actual report of BEB sites is five lines after where we see the
				# 'Bayes Empirical Bayes' string.
				next(f)
				next(f)
				next(f)
				next(f)
				next(f)
				# Process lines until we see one that starts with 'The grid'
				beb_text = next(f).strip()
				while not beb_text.startswith('The grid'):
					beb_parts = beb_text.split()
					# If we see a blank line, we will just skip it and not try to
					# parse it. We expect six parts to the BEB output line.
					if len(beb_parts) != 6:
						pass
					else:
						# Next, save any lines that have an asterisk in the third
						# field of the BEB output
						if '*' in beb_parts[2]:
							beb_sites.append(beb_parts)
					beb_text = next(f).strip()
			# When we find a line that starts with "Time used" then we have reached
			# the end of a model data output
			if line.startswith('Time used'):
				# The dN/dS value reported by PAML in the "dN & dS for each branch" section
				# is the weighted sum of w*p, which we just extracted. Treating 'w' as the
				# values to sum and 'p' as the weights for those values, we calcluate the
				# weighted mean of the various site classes.
				dn_ds = 0
				# zip() acts like 'cbind': it takes multiple lists and combines them
				# in a way where we can iterate over them together.
				# We call list(p) and list(w) because in the case of model 0, p and
				# w are single numbers, and do not work with zip() in the way that the
				# other model outputs would work. Putting them inside list() turns
				# them into lists of length 1, which works with zip()
				for weight, site_class_omega in zip(list(p), list(w)):
					site_class_contribution = float(weight) * float(site_class_omega)
					dn_ds += site_class_contribution
					# Now, stick the results into the dictionary. Key the dictionary with
					# the model name.
				if current_model in paml_dict:
					paml_dict[current_model]['w'] = w
					paml_dict[current_model]['lnL'] = lnl
					paml_dict[current_model]['np'] = np
					paml_dict[current_model]['p'] = p
					paml_dict[current_model]['dN/dS'] = str(dn_ds)
					paml_dict[current_model]['BEB.Significant'] = beb_sites
				else:
					paml_dict[current_model] = {
						'w': w,
						'lnL': lnl,
						'np': np,
						'p': p,
						'dN/dS': str(dn_ds),
						'BEB.Significant': beb_sites}
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