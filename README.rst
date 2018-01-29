Microbiome Analysis Pipeline
----------------------------

This repository includes a number of bash scripts that are used to analyze
the bacterial and fungal microbiota in human, rat, and pig skin samples, though
its utility is not limited to those microbiomes.  It is used to generate and
analyze OTU tables by processing FASTQ files from sequenced samples.

Prerequisites
-------------

This pipline was used with:

- qiime 1.9.1
- pipits 1.4.1

Additional prerequisites include:

- PEAR or VSEARCH (for joining paired-end reads)
- uclust
- usearch
- muscle (multiple alignment, ITS sequences only)

Process
-------

The directory structure for your project should look something like this:

::

  project_name
   |- sample_name1
   |   |- fastq_16S
   |   |- fastq_ITS
   |- sample_name2
   |   |- fastq_16S
   |   |- fastq_ITS
   |- sample_name3
       |- fastq_16S

The fastq_16S and fastq_ITS directories should contain the forward and
reverse reads in FASTQ format, for sequencing with 16S and ITS2 barcodes
respectively.  This pipeline assumes that barcodes have already been removed.
The fastq_16S or fastq_ITS directories should not be created if the appropriate
data is not available for the sample (as in the case of sample_name3 above).

The directories for the samples must match the sample names in the metadata
file.  All of the scripts in this repository should be run from the base
(``project_name``) directory.

**Step 1: Sample Preparation**

Running ``prepare_samples.sh`` will generate FASTA files from the FASTQ files
in the sample directories.  It uses PIPITS to reindex the reads, join the
paired-end reads (using PEAR or VSEARCH, depending on the version of PIPITS),
do quality filtering (using FASTX), and convert to FASTA.  QIIME is then used
to ensure that the sequences are named appropriately in the FASTA files.

A ``fasta`` directory will be created in each sample directory with the 16S and
ITS FASTA files.  A ``prepared_fasta_<date>`` directory will be created, which
contains a FASTA file that incorporates all 16S sequences from all samples,
and one that incorporates all ITS sequences.  All FASTA files are gzipped.
The ``scripts`` subdirectory includes the data necessary to reproduce the run,
a log file that contains the script output, and a stats file that lists the
number of sequences in the FASTA file(s) for each sample.

This script should be run as follows:

``bash prepare_samples.sh <project_name> <metadata_file_name>``

**Step 2: OTU Picking**

(script to be committed)

Uses QIIME for 16S, and PIPITS for ITS2.

- Chimeric sequences are identified and removed.
- OTU picking is run, and an OTU table is created in the biom format.
- A representative set is selected.
- Taxonomy is assigned and added to the OTU table.

**Step 3: OTU Table Analysis**

Running ``analyze_otu_table.sh`` performs analysis on the generated OTU tables
using the tools provided by QIIME.

- A phylogeny tree is created for the OTUs.
- Alpha diversity is calculated for each sample, using rarefaction, and the
  specified groups are compared.
- Beta diversity is analyzed, through jackknifing and through plots.
- Plots are generated to summarize the taxa.
- Group significance is calculated.

This script should be run as follows:

``bash analyze_otu_table.sh <project_name> <16S|ITS> <metadata_file_name> <max_rarefaction_depth> <metadata_columns>``

It should be run twice, once for 16S and once for ITS.  The metadata_columns
argument should be a comma separated list of column names from the metadata
file that define groups for comparison, such as ``TreatmentState,Description``.
To determine the value for rarefaction depth, run
``biom summarize-table -i <biom_file>``, note the Min value (the smallest
number of sequences in a sample), and choose a number somewhat smaller than
that value.

