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
- uclust
- muscle (multiple alignment)

Process
-------

(This section will be refined as the scripts are committed to the repository.)

The directory structure for you project should look something like this:

::

  project_name
   |- sample_name1
   |   |- fastq
   |- sample_name2
       |- fastq

The fastq directories should contain the forward and reverse reads in FASTQ
format, for both 16S and ITS sequencing.

#. Joining paired-end reads (script to be committed)
#. Performing OTU picking to generate OTU tables (script to be committed)
#. Run analyze_otu_table.sh to perform analysis on the generated OTU tables

