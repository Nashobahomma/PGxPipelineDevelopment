# PGx Pipeline User Manual
## Overview
## Prerequisites
### Alpine Access
Links to Alpine for Noobs slides pdf and video:
- https://olucdenver-my.sharepoint.com/personal/tonya_brunetti_cuanschutz_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Ftonya%5Fbrunetti%5Fcuanschutz%5Fedu%2FDocuments%2FActive%20HPC%20Alpine%2FAlpine%5FNoob%5FIntroduction%5Fto%5FHPC%5Fand%5FAlpine%201%2Epdf&parent=%2Fpersonal%2Ftonya%5Fbrunetti%5Fcuanschutz%5Fedu%2FDocuments%2FActive%20HPC%20Alpine&ct=1713981640547&or=OWA%2DNT%2DMail&cid=886f03e6%2D04e0%2Db856%2D85c8%2D0cfc266ce527&ga=1
- https://olucdenver-my.sharepoint.com/personal/tonya_brunetti_cuanschutz_edu/_layouts/15/stream.aspx?id=%252Fpersonal%252Ftonya_brunetti_cuanschutz_edu%252FDocuments%252FActive%20HPC%20Alpine%252Falpine_noob_04032024.mp4&ct=1713981744021&or=OWA-NT-Mail&cid=3ef9a4c1-0158-2038-ff66-d146d2aef6d4&ga=1

#### To set up an account on the computing cluster
Create an ACCESS-CI account and email rc-help@colorado.edu to request account
permissions to the ClawLab allocation. 

#### Requesting an RMACC Account on CURC Resources
1. Create an ACCESS-CI account in the ACCESS user portal: https://metrics.access-ci.org/
   - If you already have an XSEDE or ACCESS account, please do not create
   another one – just go to step 2.  
2. Email rc-help@colorado.edu, and request an account. Include the following 
information:
   - your ACCESS or XSEDE username
   - your institutional affiliation (e.g., “University of Awesome”)
   - your role (undergraduate graduate student, postdoc, staff, instructor,
     faculty or affiliated faculty)
   - your department
   - your first and last name
   - your preferred email address for communication  
3. They will set up your account and send a confirmation email.
4. Login to the RMACC OnDemand portal to access CURC resources! The first
   time you login you will be prompted to set up two-factor authentication.

#### Logging in
1. Navigate to: https://cilogon.org/authorize?client_id=cilogon:/client_id/2839765bdb5cabfdae0daa43a6614b13&scope=email%20openid%20profile%20org.cilogon.userinfo&redirect_uri=https://access-ci.org&response_type=code&state=asdfghjklkjhgfdsa
2. This page will direct you to CILogon to authenticate your Alpine session.
   - In the box "Select an Identity Provider," Select "ACCESS CI (XSEDE)"
   - Click the "Log On" box. 
3. This then takes you to the real log-on page: "ACCESS." Input your ACCESS ID
   and ACCESS password. Note: these are not the same as your Anschutz logon
   credentials.
4. This takes you to your CURC On Demand page.

#### Moving around in Alpine
1. At CURC, there is a menu bar with dropdown lists for Files, Jobs,
   Clusters, Interactive Apps, and an image for My Interactive Sessions. 
2. Files is file browser where you can view paths to data, scripts, and output 
3. Jobs is where you can keep track of the jobs that have been submitted to
   the cluster  
4. Clusters, is the dropdown list to access the Alpine shell command line.  

### Set up Conda Environment
Use the file called `PGx_Conda_Env_Package_List.txt` in the GitHub repository
to set up a Conda environment for the pipeline. Copy this file to the compute
cluster, and set up the environment in an interactive compute job:

```bash

sinteractive -N 1 -n 1 -c 4 --mem-per-cpu=4gb -t 2:00:00
module load anaconda
conda create --name CYPevol --file PGx_Conda_Env_Package_List.txt
```

**Note** the pipeline scripts are written to use an environment named
`CYPevol`. If you use a different environment name (`--name` option in the
above `conda` command), you **must** change the name in the pipeline scripts
to match.

## 1. Clone Repository (Or Download Scripts) to Cluster
1. Log in to the compute cluster with a terminal (shell).
2. Use `git` to get the latest copy of the pipeline scripts and supporting
   data files:

   ```bash
   git clone https://github.com/Nashobahomma/PGxPipelineDevelopment.git
   ```

This will create a directory named `PGxPipelineDevelopment` in your current
working directory. This directory is where all of the pipeline scripts are
stored.

## 2. Download CDS Files from NCBI
The current analysis is performed over three groups of primates. The include
6, 10, and 16 species each. The 10 species group contains the same species as
the six species group with four species added and the 16 species group
likewise but with six species added on top of the 10 species group. Each
larger group contains the addition of species evolutionarily more distant
from humans. That is to say that the six species group is species closes to
humans and the 16 species group includes species as far out as lorises and
tarsiers.

This analysis pipeline uses only CDS (coding sequences) fasta files.

CDS files of each species' Reference Sequence need to be updated and thus
analysis rerun once per quarter (4x per year) as NCBI updates (adds and/or
removes) records.

Sequences need to first be downloaded to your local desktop before
extracting the necessary fasta file for upload to the google cloud server.

Use Chrome or Mozilla to download sequences. Safari use will present
problems: https://www.ncbi.nlm.nih.gov/datasets/docs/reference-docs/mac-zip-bug/

When downloading new CDS files from NCBI:

1. Go to https://www.ncbi.nlm.nih.gov/
2. In the search window, enter the binomial scientific name of a species
   (e.g. Pan trodlodytes).
3. Click on the link with the name of that species. In the "Genome," then 
   "Reference Genome" section, select the "Download" button.
4. In the "Download Package" pop-up window:
    - Click on the "RefSeq only" radio dial
    - Uncheck "Genomic sequences, (FASTA)" and check “Genomic Coding
      Sequences (FASTA)”
    - Rename the GCF zip file at the bottom of the pop-up to reflect the
      scientific name of the species and the name of the assembly (should
      be the reference assembl) (e.g. Mmulatta_Mmul_10)
5. Finally, click download, and place the file in a local (your computer)
    directory that you create called, CDS_FASTA
6. Unzip each file. The fasta file will be in ncbi_dataset/data/GCF_XXXXX;
   it will have an .fna extension.
7. Save the .fna file, change the extension name to .fasta (e.g.
    Pan_trodlodytes.fasta, note somewhere the reference sequence name and
    date of download).
8. If the CDS sequences downloaded do not have a peptide product annotated
    for a gene (sequence), it will not be included in the CDS file.

## 3. Edit `Lemma.sh` (Configure Pipeline)
The file located at `Final_Pipeline_Scripts/Lemma.sh` relative to the cloned
GitHub repository is file file that you must edit to configure your pipeline
run. Use the Open OnDemand file edtior or a command line text editor to make
changes to match your analysis. These variables should be edited:

- `_PIPE_SCRIPTS_FROM_GITHUB`: Set to the full path of the location on the cluster
  that the GitHub repository was cloned into.
- `_PIPE_EMAIL_ADDRESS`: Set to the email address to where you will receive
  e-mail updates.
- `_PIPE_ALL_DATA`: Set to the path where "static reference data" (e.g., NCBI
  CDS files) and supporting data files, including those generated during the
  pipeline run. In practice, it is set to the *parent* directory of the NCBI
  CDS files directory.
- `_PIPE_COHORT_MEMBERS`: Set to an array of FASTA files that correspond to
  the species that should be part of the analysis. It is flexible in that
  species can be excluded by "commenting" (prefixing with `#`) the corresponding
  FASTA file name. Make sure it is defined with only one file name per line,
  and that the parentheses `()` are not commented out.

### 3.1: Optional Configuration Parameters
These values are derived from other user-configurable values in `Lemma.sh`,
but you can override them. It is generally do not recommended unless you want
to have a highly customized output structure.

- `_PIPE_RUN_NICKNAME`: Set to name of the desired output folder. By default,
  it is the number of species and the execute date, separated by an underscore.
- Job resource request parameters (e.g., partition, walltime, cores, memory)

## 4. Run `Palea.sh` (Execute Pipeline)
Navigate to the `PGxPipelineDevelopment/Final_Pipeline_Scripts` directory. Run
the `Palea.sh` file with `bash`:

```bash
cd /projects/anashoba@xsede.org/PipelineFiles/PGxPipelineDevelopment/Final_Pipeline_Scripts
bash Palea.sh
```

**Note** that the `Palea.sh` script depends on all of the scripts from the
GitHub repository being in the same directory. Do not move the scripts around
before trying to run the pipeline.

### 4.1 Pipeline Output Files
The output is written to the directory specified by `_PIPE_ALL_DATA`, in a
sub-directory named according to the value of `_PIPE_RUN_NICKNAME`. The
following files and directories are produced:

- `Checkpoints`: Directory of checkpoint files for tracking pipeline progress
- `PGx_Pipeline_Execution_Record.txt`: Text file wtih details of who ran the
  pipeline, when, from which directory, and the job IDs of the pipeline jobs.
- `Scheduler_Logs`: Directory of slurm scheduler log files
- `Species_list.txt`: Text file that lists which species were analyzed
- `Step_00_Orthofinder_TargetOGs`: Orthofinder groups that have human CYPs of
  interest
- `Step_01_PAML_Seq_Inputs`: Aligned FASTA files for PAML input
- `Step_02_PAML_Gene_Trees`: Gene trees for PAML input
- `Step_03_PAML_Control_Files`: Control files that specify PAML models for each
  gene group
- `Step_04Acc_PAML_Accessory_Files`: "Extra" files produced by PAML during runs
- `Step_04_PAML_Runs`: Main output files from PAML (omegas and lnL values)

## 5. Analyze Results
The output table produced has these four fixed columns:
- `Orthogroup.ID`: Orthofinder ID for the ortholog group
- `CYP.Name`: Gene symbol of the human CYP of interest
- `Human.Protein.ID`: NCBI protein ID for the human CYP of interest
- `Num.Species`: Number of species in the alignment analyzed by PAML

Additionally, it produces these columns for each **model** run by PAML:
- `lnL`: :og-likelihood (natural log) of the model
- `dNdS`: **Gene-wide average** dN/dS ratio
- `np`: Number of parameters in the PAML model (necessary for chi-squared test)

Neutral model (0,1,7,8a) additionally have these columns:
- `W`: Omega values for the site categories indentified by PAML
- `P`: Proportion of the gene in each site category

**Note**: models that fit multiple site categories will have each category
value joined by semicolons (`;`).

Selection models (2, 8) have these columns:
BEB: Bayes-empirical-Bayes (cite here): these are amino acids that have
statistically significant evidence of positive selection (w>1, P<0.05 [Note,
different from `P` in the neutral model output])

- `BEB.Position`: Significant BEB alignment position
- `BEB.AminoAcid`: Significant BEB amino acid
- `BEB.Probability`: P-value for BEB
- `BEB.W`: Omega estimate for the significant position
- `BEB.WSE`: Standard error of the estimated omega

**Note**: when multiple sites cross the threshold for significance, they are
joined by a semicolon (`;`) in the output CSV.

### 5.1: Likelihood Ratio Tests for Model Selection
Perform the following comparisions for model selection:

- Model 2 (selection) vs Model 1 (neutral)
- Model 8 (selection) vs Model 7 (neutral)
- Model 8 (selection) vs Model 8a (netural)

To calculate the likelihood ratio test score and P-value, use the Chi-squared
distribution. In R, for example model 2 vs model 1:

```r
test_stat <- 2*(Model2.lnL - Model1.lnL)
test_pval <- pchisq(test_stat, df=Model2.np-Model1.np, lower.tail=FALSE)
```

If `test_pval` is less than 0.05 (or your favorite threshold), then you can
reject the netural model (null model).

Perform these comparisons for each orthologous gene group, and choose the model
that has best statistical support.