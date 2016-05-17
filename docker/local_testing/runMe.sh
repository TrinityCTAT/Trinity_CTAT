#!/bin/bash

set -e

export DISCASM_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/DISCASM
export STAR_FUSION_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/STAR-Fusion
export FUSION_INSPECTOR_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/FusionInspector
export TRINITY_HOME=/home/unix/bhaas/GITHUB/trinityrnaseq

utildir=`dirname $0`/..

export PERL5LIB=${utildir}/PerlLib:${PERL5LIB}


CTAT_LIB=/seq/RNASEQ/public_ftp/CTAT_lib.tar.gz
LEFT_FQ=/seq/RNASEQ/CANCER_FUSION_RESOURCES/BT-474/BT474.Left.fq.gz
RIGHT_FQ=/seq/RNASEQ/CANCER_FUSION_RESOURCES/BT-474/BT474.Right.fq.gz

cmd="${utildir}/util/CTAT_fusion_wrapper.pl $LEFT_FQ $RIGHT_FQ $CTAT_LIB DISCASM FusionInspector"

eval $cmd

