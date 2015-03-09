#!/bin/bash -l

umask 0002

. /broad/software/scripts/useuse

reuse GCC-4.4
reuse LSF
reuse Perl-5.8
reuse Python-2.7
reuse  Java-1.7
reuse .coreutils-8.22
reuse .samtools-0.1.19

export PYTHONPATH=/seq/regev_genome_portal/SOFTWARE/PYTHON_LIB:${PYTHONPATH}
export HISAT_HOME=/seq/regev_genome_portal/SOFTWARE/HISAT/current
export TRINITY_HOME=/home/unix/bhaas/GITHUB/trinityrnaseq
export PATH=/home/unix/bhaas/utilities:${PATH}

./runMe.pl



