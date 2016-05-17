#!/bin/bash

set -e

export DISCASM_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/DISCASM
export STAR_FUSION_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/STAR-Fusion
export FUSION_INSPECTOR_HOME=/home/unix/bhaas/GITHUB/CTAT_FUSIONS/FusionInspector
export TRINITY_HOME=/home/unix/bhaas/GITHUB/trinityrnaseq

utildir=`dirname $0`/..

export PERL5LIB=${utildir}/PerlLib:${PERL5LIB}

cmd="${utildir}/util/CTAT_fusion_wrapper.pl $STAR_FUSION_HOME/testing/reads_1.fq.gz \
                                            $STAR_FUSION_HOME/testing/reads_2.fq.gz \
                                           /seq/RNASEQ/public_ftp/CTAT_lib.tar.gz \
                                           DISCASM FusionInspector"

eval $cmd

