# Tasks for Pipeline Development
Date: 2024-05-19

## Need
- will have to put into the trans-back-ogsel .sh the location of the orthofinder generated .tsv file --- as well as properly referencing the main versions of the orthofinder output .fa files
- make sure that the directories for each Step_XX are generated in the scripts  somewhere so that the end user does not need to generate these manually
- to address the working environment.. currently using one cloned from kendra ferrier... tom has stragedies for this
- in the write ctl file .py script need to generalize the output directory on line 34

## Want
- for the sbatch directive pertaining to speficying date run of .err and .out files ==> auto insert the date that overcoat script is run
- in the make gene trees step... specify a stub of sorts indicating at what point/where/etc the foreground branch needs to be labeled for the use of this pipeline to run PAML codeml branch-site models
- R script to compare the likelihoods from the PAML models and to parse out the relevant omega values
