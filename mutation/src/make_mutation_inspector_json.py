#!/usr/bin/env python

import argparse
import csv
import json
import os

prsr_arguments = argparse.ArgumentParser( prog = "make_mutation_inspector_json.py", description = "Creates a json object for the mutation pipeline to be displayed in the mutation inspector app.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--sample", dest="str_sample_name", required=True, help = "Name to display for sample or run." )
prsr_arguments.add_argument( "--tab", dest="str_cancer_tab", required=True, help = "Path to cancer tab file (Holds variants of interest)." )
prsr_arguments.add_argument( "--bam", dest="str_bam", required=True, help = "Path to bam that variants were called from." )
prsr_arguments.add_argument( "--bam_index", dest="str_bai", required=True, help = "Path to bai for bam that variants were called from." )
prsr_arguments.add_argument( "--bed", dest="str_bed", required=True, help = "Path to bed that annotate the reference genome." )
prsr_arguments.add_argument( "--bed_index", dest="str_bed_index", required=True, help = "Path to the index file (eg. idx) for the bed file." )
prsr_arguments.add_argument( dest="str_output_json_file", help = "Output json file." )
args = prsr_arguments.parse_args()

# Constants
C_STR_BAM = "BAM"
C_STR_BAM_INDEX = "BAM_INDEX"
C_STR_BED = "BED"
C_STR_BED_INDEX = "BED_INDEX"
C_STR_SAMPLE = "SAMPLE"
C_STR_SNV = "SNV"
C_STR_DELIMITER = "\t"

dict_json = { C_STR_BAM : args.str_bam,
              C_STR_BAM_INDEX : args.str_bai,
              C_STR_BED : args.str_bed,
              C_STR_BED_INDEX : args.str_bed_index,
              C_STR_SAMPLE : os.path.basename( args.str_sample_name ) }
ldict_entries = []

lstr_header = None

# Read in the cancer tab file
with open( args.str_cancer_tab, "r" ) as hndl_in:
  for lstr_line in csv.reader( hndl_in, delimiter = C_STR_DELIMITER ):
    if lstr_header is None:
      lstr_header = lstr_line
      continue
    ldict_entries.append( dict( zip( lstr_header, lstr_line ) ) )
dict_json[ C_STR_SNV ] = ldict_entries

# Write pretty json to file
with open( args.str_output_json_file, "w" ) as hndl_out:
  hndl_out.write( json.dumps( dict_json, sort_keys=True, indent=2 ) )
