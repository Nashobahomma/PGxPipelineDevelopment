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


## 2024-06-02

Finished a first draft of all pipeline steps. Currently testing it on Alpine. Copy of the test execution record:

```
Record of PGx CPY450 PAML pipeline.
Pipeline version: 0.0.0
Execution time: 2024-06-02_14:55:53
Executing user: anashoba@xsede.org
Nickname of run: 16_2024-06-02
Working directory: /projects/anashoba@xsede.org/PipelineFiles/PGxPipelineDevelopment/Final_Pipeline_Scripts
Invoking command: Palea.sh  ''
Input CDS FASTA directory: /projects/anashoba@xsede.org/PipelineFiles/NCBI_CDS_Files
Number of CDS FASTA files: 16
Final output directory: /projects/anashoba@xsede.org/PipelineFiles/16_2024-06-02
Step 00: Run_Orthofinder has job ID 6230121
Step 01: Prepare_PAML_Sequences has job ID 6230122
Step 02: Make_Gene_Trees has job ID 6230123
Step 03: Make_PAML_Control_Files has job ID 6230124
Step 04: Run_PAML_Site_Models has job ID 6230125
```

Next steps:

- Smooth out wrinkles in pipeline execution.
- Write user guide
- Amber work on slides for 17 June presentation: will add data as it becomes available
also next time, run the six and ten species versions
Add data - smooth slides
