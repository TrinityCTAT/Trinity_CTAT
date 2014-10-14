Requirements
============

The following software is required before running the RNAseq variant calling pipeline:

* GATK: http://www.broadinstitute.org/gatk
* Java 1.7
* PICARD tools: http://picard.sourceforge.net/command-line-overview.shtml
* Python >= 2.7
* R: http://cran.r-project.org/
* SciEDPiper: https://github.com/SciEDPipeR/SciEDPipeR
* STAR Aligner 2.3.0: http://code.google.com/p/rna-star

Current versions
=================
At the time of writting this documentation the following versions of software
are being used successfully. Other versions of the software dependencies may
also be viable.

* GATK
** GenomeAnalysisTK-3.1-1-g07a4bf8
* Java
** java version "1.7.0_51"
** Java(TM) SE Runtime Environment (build 1.7.0_51-b13)
** Java HotSpot(TM) 64-Bit Server VM (build 24.51-b03, mixed mode)
* PICARD tools
** AddOrReplaceReadGroups.jar
*** Version: 1.764(96ec474b689a429463e04256383babd1f62efd88_1410735622)
** MarkDUplicates.jar
*** Version: 1.764(96ec474b689a429463e04256383babd1f62efd88_1410735622)
* Python
** Python 2.7.1 (r271:86832, Apr 17 2012, 22:46:32) 
** [GCC 4.4.3] on linux2
* R
** R version 3.0.2 (2013-09-25) -- "Frisbee Sailing" (64-bit)
* SciEDPiper
** v0.1.0
* STAR Aligner
** 2.3.0


Make sure to append to your PATH the path the STAR aligner and SciEDPiper.
Also, make sure Picard-Tools and GATK jars are available.


Quick Start
===========

You will need the following files:
* A reference genome ( fasta )
* A vcf file of all know variants for your reference genome ( vcf, for example a dbSNP vcf )
* 2 paired RNASeq samples ( fasta or fasta.qz)

Either put the rnaseq_mutation_pipeline in your path or go to its directory.

Use the following command:
Note: files in {} will need to be updated with actual files.

python rnaseq_mutation_pipeline.py --reference {reference_genome.fa} --vcf {reference_genome.vcf} --left {left_sample.fa} --right {right_sample.fa} --threads 8 --log log.txt --out_dir rnaseq_mut_out


Pipeline behavior
=================

Note that if an error occurs on a command, the pipeline should not produce files which would have resulted 
from that command, so no partial files should be created in the pipeline. As well, if the pipeline was stopped
midway for some reason, when ran again with the same settings, the pipeline will pick-up where it left off.

There is an option to run the pipeline with a clean command line argument. This is not a default behavior.
If used, when an intermediary file is made in the pipeline which is not an input file and not a final result,
it is only kept until it is not longer needed for commands and then it is deleted. 
This keeps the foot print of the run down.

If any of these behaviors are not true, please let me know ttickle@broadinstitute.org


LSF Platform integration
========================

Each step of the pipeline can be sequentially executed using bsub. Set the Bsub queue for the job using -b/bsub_queue. Optionally, memory
requirements can be given for LSF node selection using -m/max_bsub_memory .


Gzipped files
=============

The .gz extension is recognized by the pipeline and any commands using this will automatically have the .gz file unzipped and piped in for
use making any command compatible with zipped files. Please use them space footprints.


Updating paths
==============

Several jars are used in this pipeline. You may supply the pipeline to update any command. This occurs with the -u/update_command .


Reducing processing time
========================

The indexing step takes is the longest step and can create large files ( GBs ). This happens twice, the 
first indexing is a only dependent on the reference genome and so can be shared in the processing of many samples.
It is possible to index the reference genome first and then supply the directory of the reference genome index
to the pipeline.

We have made a script to make building the index easier.
Here is an example command on building the index.
Note: Files in {} need to be updated with actual files

python rnaseq_mutation_indexer.py --reference {reference_genome.fa} --threads 8 --log indexing_log.txt --out_dir premade_index_directory

If you have an index already made you can then supply the rnaseq_mutation_pipeline.py with the location 
of the index using the optional argument -i / --index

python rnaseq_mutation_pipeline.py --reference {reference_genome.fa} --vcf {reference_genome.vcf} --left {left_sample.fa} --right {right_sample.fa} --threads 8 --log log.txt --out_dir rnaseq_mut_out --index {premade_index_directory}


Argument Help
=============

More specific help with command line arguments can be found with the following
command:

python rnaseq_mutation_indexer.py --help
