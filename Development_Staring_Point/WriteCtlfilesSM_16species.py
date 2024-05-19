 #!/usr/bin/env python

"""REWORKED 10Oct2022 for Site models ...(1,2,7,8,8a?)Code to Generate PAML control file. Adapted from code 
by KClaw by A.Nashoba and T.Kono August 2022. This code will generate
a control file for one alignment. A separate script will be
written to batch this action across multiple alignments. 

Usage:
 
python Write_Ctl_files_Site_Models.py /path/to/OGXXXXX_RepOrthologues.fa /path/to/OGXXXXX.raxml.bestTree NSSITES > /path/to/OGXXXXX_ctl.txt

Where NSSITES is '8a' or '01278', to write control files for the null model (8a) or the
alternate model, respectively.
"""

import sys
import os

try:
	og_in = sys.argv[1]
	tree_in = sys.argv[2]
	nssites_type = sys.argv[3]
except IndexError:	
	sys.stderr.write(__doc__ + '\n')
	sys.exit(1)


# DEFINE CONSTANTS FOR THE PAML CONTROL FILE
# Edit these to adjust the control file that gets written for each orthogroup
# We are defining these parameters as constants because they will be static for
# each control file that gets written. 
#	OUTPUT PATH # !This needs to be generated outside a script (e.g. in Finder) before running PAML.
#PAML_OUTPUT_DIRECTORY = '/Users/ambereule-nashoba/Desktop/Dropbox/KClaw/PipelineB/16species_57Genes/16species_57Genes_Workflow/Step_06_PAML_Runs'
PAML_OUTPUT_DIRECTORY='/projects/anashoba@xsede.org/ForAlpine/Step_06_PAML_Runs'
#	SCREEN OUTPUT SETTINGS
NOISY = '9'
VERBOSE = '0'
#	GAMMA DISTRIBUTION FOR SUBSTITUTION RATE AMONG SITES
#		-3 = Pairwise Bayesian
#		-2 = ML by numerical iteration
#		0 = discrete gamma
# Runmode = 0: user-specified tree file (i.e., do not estimate a new tree from alignment)
RUNMODE = '0'
#	INPUT SEQUENCE TYPE
SEQTYPE = '1'
CODONFREQ = '1'
#	AMINO ACID DISTANCE MATRIX
AADIST = '0'
#	MODEL CLASS
#	See page 35 in pamlDOC.pdf
#	These will not be specified as constants because they will change depending
# 	on user input.
#	Null for branch-site model is model=2; NSsites=2; fix_omega=1; omega=1
#	GENETIC CODE
ICODE = '0'
#	SEQUENCE PARTITIONS
MGENE = '0'
#	KAPPA
FIX_KAPPA = '0'
KAPPA = '2'
#	ALPHA
FIX_ALPHA = '1'
ALPHA = '0'
MALPHA = '0'
#	NUMBER OF SITE CATEGORIES 
NCATG = '10'
#	CLOCK
CLOCK = '0'
#	STANDARD ERROR
GETSE = '0'
#	RATE ANCESTOR
RATEANCESTOR = '1'
#	SMALL DIFFERENCE
SMALL_DIFF = '.5e-6'
#	CLEAN DATA
#		0 = Include alignment columns with '?' characters in the analysis
#		1 = Exclude alignment columns with '?' characters in the analysis
CLEAN_DATA = '0'
#	NUMBER OF DATASETS PER INPUT FILE
NDATA = '1'
#	ESTIMATION METHOD
#		0 = Estimate all parameters simultaenously
#		1 = Estimate parameters one at a time
METHOD = '0'
# For site models, we set MODEL=0
MODEL = '0'
# The second argument will specify the NSsites value. 
# Now, use the second argument (NULL/ALT) to set the PAML control options
# for specifying the null or the alternate models to fit.
if nssites_type == '01278':
	# For NSsites 01278, we set NSsites=0 1 2 7 8; fix_omega=0; omega=0
	NSSITES = '0 1 2 7 8'
	FIX_OMEGA = '0'
	OMEGA = '0'
elif nssites_type == '8a':
	# For 8a (a null hypothesis), we set NSsites=2; fix_omega=0
	NSSITES = '8'
	FIX_OMEGA = '1'
	OMEGA = '1'
else:
	sys.stderr.write('Error! NSSITES must be one of "8a" or "01278" (case sensitive).\n')
	sys.stderr.write('The NSSITES you specified was ' + nssites_type + '\n')
	sys.exit(1)

# Next, use the input filename to designate the input and output filenames
SEQFILE = os.path.abspath(os.path.expanduser(og_in))
TREEFILE = os.path.abspath(os.path.expanduser(tree_in))
#	Isolating the orthogroup name from the input alignment filename
og_name = os.path.basename(og_in).split('_')[0]
# 	Building the output path from the PAML_OUTPUT_DIRECTORY constant and the OG name
#	Include the hypothesis/model type (NULL/ALT) in the output name
OUTFILE = os.path.join(PAML_OUTPUT_DIRECTORY, og_name + '.' + nssites_type + '.out')

# Finally, write the control file to standard output channel
print('seqfile = ' + SEQFILE)
print('treefile = ' + TREEFILE)
print('outfile = ' + OUTFILE)
print('noisy = ' + NOISY)
print('verbose = ' + VERBOSE)
print('runmode = ' + RUNMODE)
print('seqtype = ' + SEQTYPE)
print('CodonFreq = ' + CODONFREQ)
print('aaDist = ' + AADIST)
print('model = ' + MODEL)
print('NSsites = ' + NSSITES)
print('icode = ' + ICODE)
print('Mgene = ' + MGENE)
print('fix_kappa = ' + FIX_KAPPA)
print('kappa = ' + KAPPA)
print('fix_omega = ' + FIX_OMEGA)
if nssites_type == '8a':
	print('omega = ' + OMEGA)
print('fix_alpha = ' + FIX_ALPHA)
print('alpha = ' + ALPHA)
print('Malpha = ' + MALPHA)
print('ncatG = ' + NCATG)
print('clock = ' + CLOCK)
print('getSE = ' + GETSE)
print('RateAncestor = ' + RATEANCESTOR)
print('Small_Diff = ' + SMALL_DIFF)
print('cleandata = ' + CLEAN_DATA)
print('ndata = ' + NDATA)
print('method = ' + METHOD)
