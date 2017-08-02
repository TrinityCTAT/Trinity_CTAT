#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import subprocess
import argparse
import logging

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="firecloud fusion pipeline runner", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--left_fq", type=str, required=True, help="left fq file")

    parser.add_argument("--right_fq", type=str, required=True, help="right fq file")

    parser.add_argument("--genome_lib_tar_gz", required=True, help="ctat genome lib as tar.gz")

    parser.add_argument("--output_dir_name", required=True, help="output name")

    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")

    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)



    # prep the genome lib
    subprocess.check_call("tar xvf {}".format(args.genome_lib_tar_gz), shell=True)

    genome_lib_dir = re.sub(".tar.gz$", "", args.genome_lib_tar_gz)


    ## run star-fusion

    outdir = args.output_dir_name
    
    starF_cmd = str("/usr/local/src/STAR-Fusion-v1.0.0/STAR-Fusion" +
                    " --left_fq {}".format(args.left_fq) +
                    " --right_fq {}".format(args.right_fq) +
                    " --genome_lib_dir {}".format(args.genome_lib_dir) +
                    " --extract_fusion_reads " +
                    " --FusionInspector validate " +
                    " --denovo_reconstruct " +
                    " --annotate " +
                    " --examine_coding_effect " +
                    " --output_dir {}".format(args.output_dir_name)
                    )

    subprocess.check_call(starF_cmd, shell=True)
    
    ## package up outputs
    cmd = str("tar -cvfz ${}.tar.gz".format(outdir) +
              " {}/star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv".format(outdir) +
              " {}/star-fusion.fusion_evidence_reads_1.fq ".format(outdir) +
              " {}/star-fusion.fusion_evidence_reads_2.fq ".format(outdir) +
              " {}/FusionInspector/finspector.spanning_reads.bam ".format(outdir) +
              " {}/FusionInspector/finspector.junction_reads.bam ".format(outdir) +
              " {}/FusionInspector/finspector.fusion_predictions.final.abridged.FFPM ".format(outdir) +
              " {}/FusionInspector/finspector.gmap_trinity_GG.fusions.gff3 ".format(outdir) +
              " {}/FusionInspector/finspector.gmap_trinity_GG.fusions.fasta ".format(outdir)
              )

    subprocess.check_call(cmd, shell=True)


    # done.
    sys.exit(0)

####################

if __name__ == "__main__":
    main()
    


