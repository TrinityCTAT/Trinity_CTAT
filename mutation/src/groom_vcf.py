#!/usr/bin/env python

# Constants
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"

import argparse
import csv

prsr_arguments = argparse.ArgumentParser( prog = "groom_vcf.py", description = "Cleans VCF files for bcftools.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "str_input_file", help = "Input vcf file." )
prsr_arguments.add_argument( "str_output_file", help = "Output groomed vcf file." )
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []

# Read in vcf file
if args.str_input_file:
  with open( args.str_output_file, "w" ) as hndl_out:
    with open( args.str_input_file, "r" ) as hndl_vcf:
      for lstr_line in csv.reader( hndl_vcf, delimiter = STR_VCF_DELIMITER ):
      
        # Work with comments
        if lstr_line[0][0] == CHR_COMMENT:
          # Make sure there are not spaces (especially at the end).
          lstr_line = [ str_token for str_token in lstr_line if str_token ]

          # Store
          lstr_vcf.append( STR_VCF_DELIMITER.join( lstr_line ) )
          continue

        # Work with the body of the file

        # Store
        lstr_vcf.append( STR_VCF_DELIMITER.join( lstr_line ) )

    for str_out_line in lstr_vcf:
      hndl_out.write( str_out_line + "\n" )
