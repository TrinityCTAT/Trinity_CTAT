#!/usr/bin/env python

import argparse
import csv
import json

prsr_arguments = argparse.ArgumentParser( prog = "make_mutation_inspector_json.py", description = "Creates a json object for the mutation pipeline to be displayed in the mutation inspector app.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "str_cancer_tab", help = "Path to cancer tab file (Holds variants of interest)." )
prsr_arguments.add_argument( "str_bam", help = "Path to bam that variants were called from." )
prsr_arguments.add_argument( "str_bai", help = "Path to bai for bam that variants were called from." )
prsr_arguments.add_argument( "str_output_json_file", help = "Output json file." )
args = prsr_arguments.parse_args()

# Constants
C_STR_BAM = "BAM"
C_STR_BAM_INDEX = "BAM_INDEX"
C_STR_SAMPLE = "SAMPLE"
C_STR_SNV = "SNV"

dict_json = { C_STR_BAM : args.str_bam,
              C_STR_BAM_INDEX : args.str_bai,
              C_STR_SAMPLE : os.path.splitext( os.path.basename(args.str_bam) )[0] }
ldict_entries = []

lstr_header = None

# Read in the cancer tab file
with open( args.str_cancer_tab, "r" ) as hndl_in:
  for lstr_line in csv.reader( hndl_in, delimiter = STR_TAB_DELIMITER ):
    if lstr_header is None:
      lstr_header = lstr_line
      continue
    ldict_entries.append( dict( zip( lstr_header, lstr_line ) ) )
dict_json[ C_STR_SNV ] = ldict_entries

# Write pretty json to file
with open( args.str_output_json_file, "w" ) as hndl_out:
  hndl_out.out( json.dumps( dict_json, sort_keys=True, indent=2 ) )
