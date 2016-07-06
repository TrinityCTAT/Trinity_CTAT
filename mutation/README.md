# Overview

Mutation detection in RNA-Seq highlights the GATK Best Practices in RNA-Seq variant calling, several sources of variant annotation, and filtering based on CRAVAT.

# Public and Free Galaxy Instance

This tool, as well as others, is installed at Indiana University on a publically available Galaxy instance free to use for cancer research. Please register and use this service (https://galaxy.ncgas-trinity.indiana.edu/user/login). To run the tool locally continue reading.

# Quick Start

Download the resource files found at https://data.broadinstitute.org/Trinity/CTAT/mutation/mutation_resources.tar.gz

Download demo data at https://data.broadinstitute.org/Trinity/CTAT/mutation/demo_data.tar.gz

You will need the following files ( found in the resource bundle ):
* A reference genome ( fasta )
* A vcf file of all know variants for your reference genome ( vcf, for example a dbSNP vcf )
* 2 paired RNASeq samples ( fasta or fasta.qz)

Either put the rnaseq_mutation_pipeline in your path or go to its directory.

Use the following command:
Note: values in {} will need to be updated with actual files / options.

```
python rnaseq_mutation_pipeline.py \
    --reference {reference_genome.fa} \
    --vcf {reference_genome.vcf} \
    --left {left_sample.fa} \
    --right {right_sample.fa} \
    --threads 8 \
    --cosmic_vcf {cosmic.vcf} \
    --darned {darned.txt} \
    --radar {radar.txt} \
    --log log.txt \
    --tissue_type {tissue} \
    --email {your@email.com} \
    --cravat_annotation_header { Trinity_CTAT/mutation/headers/cravat_annotation.txt } \
    --is_hg19 \
    --out_dir rnaseq_mut_out \
    --plot \
    --index {star_index} \
    --bed {reference.bed}
```

Using the demo data and the included demo resource files and assuming all files are in your current directory the command would be:
```
python rnaseq_mutation_pipeline.py \
    --reference Hg19_11.fa \
    --vcf dbsnp_FLI1.vcf \
    --left FLI1.left.fq \
    --right FLI1.right.fq \
    --threads 8 \
    --cosmic_vcf CosmicCodingMuts_v72_grch37_updated_chr.vcf.gz \
    --darned darned_hg19.txt \
    --radar radar_hg19_v2.tx
    --log log.txt \
    --tissue_type Blood-Lymphocyte \
    --email change@email.com \
    --cravat_annotation_header cravat_annotation.txt \
    --is_hg19 \
    --out_dir rnaseq_mut_out \
    --plot \
    --index Hg19_11 \
    --bed hg19.refGene.sort.bed
```

The following arguments are optional and needed to turn on CRAVAT associated functionality:
* --tissue_type
* --email
* --cravat_annotation_header
* --is_hg19

The following arguments are optional and needed for RNA editing filtering:
* --darned
* --radar

The following arguments are needed for some of the functionality involved in annotation and filtering:
* --cosmic_vcf
* (Please download from http://cancer.sanger.ac.uk/cosmic, bgzip, and make sure Chr nomenclature matches).
* --vcf

# Requirements

The following software is required before running the RNA-Seq variant calling pipeline. Please look below for the supported versions. PLEASE use the supported versions of the tools, other versions are not supported and may not work:

# Current versions

At the time of writing this documentation the following versions of software
are being used successfully. Other versions of the software dependencies may
also be viable and are not supported.

* GATK
  * GenomeAnalysisTK-3.1-1-g07a4bf8
  * http://www.broadinstitute.org/gatk
* Java
  * java version "1.7.0_51"
  * Java(TM) SE Runtime Environment (build 1.7.0_51-b13)
  * Java HotSpot(TM) 64-Bit Server VM (build 24.51-b03, mixed mode)
* PICARD tools
  * AddOrReplaceReadGroups.jar
    * Version: 1.764(96ec474b689a429463e04256383babd1f62efd88_1410735622)
  * MarkDUplicates.jar
    * Version: 1.764(96ec474b689a429463e04256383babd1f62efd88_1410735622)
  * http://broadinstitute.github.io/picard
* Python
  * Python 2.7.1 (r271:86832, Apr 17 2012, 22:46:32) 
  * [GCC 4.4.3] on linux2
* R
  * R version 3.0.2 (2013-09-25) -- "Frisbee Sailing" (64-bit)
  * http://cran.r-project.org
* SciEDPiper
  * Commit a96c6ec199c4c97f86c4a9979bd18baa1759ae91
  * https://github.com/SciEDPipeR/SciEDPipeR
* STAR Aligner
  * 2.3.0
  * https://github.com/alexdobin/STAR

Make sure to append to your PATH the path the STAR aligner and SciEDPiper.
Also, make sure Picard-Tools and GATK jars are available.


# Resources:

Download the resource files for HG19 at
https://data.broadinstitute.org/Trinity/CTAT/mutation/mutation_resources.tar.gz


# Pipeline behavior

Note that if an error occurs on a command, the pipeline should not produce files that would have resulted 
from that command, so no partial files should be created in the pipeline. As well, if the pipeline was stopped
midway for some reason, when ran again with the same settings, the pipeline will pick-up where it left off.

There is an option to run the pipeline with a clean command line argument. This is not a default behavior.
If used, when an intermediary file is made in the pipeline that is not an input file and not a final result,
it is only kept until it is no longer needed for commands and then it is deleted. 
This keeps the footprint of the run down.

If any of these behaviors are not true, please let me know ttickle@broadinstitute.org


# Gzipped files

The pipeline recognizes the .gz extension and any commands using this will automatically have the .gz file unzipped and piped in for
use making any command compatible with zipped files. Please use them space footprints.


# Updating paths

Several jars are used in this pipeline. If you need, you may supply the path of any command. This occurs with the -u/update command .


# Argument Help

More specific help with command line arguments can be found with the following
commands:

python rnaseq_mutation_pipeline.py --help
