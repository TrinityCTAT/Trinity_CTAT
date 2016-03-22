#!/bin/bash

set -e

export DISCASM_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/DISCASM
export STAR_FUSION_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/STAR-Fusion
export FUSION_INSPECTOR_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/FusionInspector

utildir=`dirname $0`/..

export PERL5LIB=${utildir}/PerlLib:${PERL5LIB}

cmd="${utildir}/util/CTAT_fusion_wrapper.pl /seq/RNASEQ/CANCER_FUSION_RESOURCES/BT-474/BT474.Left.fq.gz \
                                           /seq/RNASEQ/CANCER_FUSION_RESOURCES/BT-474/BT474.Right.fq.gz \
                                           /seq/RNASEQ/public_ftp/CTAT_lib.tar.gz \
                                           DISCASM"

eval $cmd

