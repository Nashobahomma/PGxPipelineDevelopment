# PGx Pipeline User Manual
## Overview
## Prerequisites
### Alpine Access
## 1. Clone Repository (Or Download Scripts) to Cluster
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
## 4. Run `Palea.sh` (Execute Pipeline)
## 5. Analyze Results