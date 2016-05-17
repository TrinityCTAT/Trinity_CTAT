#!/bin/bash

set -e

export DISCASM_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/DISCASM
export STAR_FUSION_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/STAR-Fusion
export FUSION_INSPECTOR_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/FusionInspector
export TRINITY_HOME=/home/unix/bhaas/GITHUB/trinityrnaseq

utildir=`dirname $0`/..

export PERL5LIB=${utildir}/PerlLib:${PERL5LIB}


CTAT_LIB=/seq/RNASEQ/public_ftp/GRCh37_gencode_v19_FULL.tgz
LEFT_FQ=$STAR_FUSION_HOME/testing/reads_1.fq.gz
RIGHT_FQ=$STAR_FUSION_HOME/testing/reads_2.fq.gz

cmd="${utildir}/util/CTAT_fusion_wrapper.pl $LEFT_FQ $RIGHT_FQ $CTAT_LIB DISCASM FusionInspector"

eval $cmd

