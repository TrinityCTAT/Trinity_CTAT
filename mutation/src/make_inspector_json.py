#!/usr/bin/env python

import argparse
import csv
from operator import itemgetter
import json
import os

# Constants
c_STR_INSPECTOR_TP = "TP"
c_STR_INSPECTOR_FP = "FP"
c_STR_INSPECTOR_FN = "FN"
c_STR_INSPECTOR_RNA_BAM = "RNA"
c_STR_INSPECTOR_DNA_BAM = "DNA"
c_STR_INSPECTOR_RNA_VCF = "RNA_VCF"
c_STR_INSPECTOR_DNA_VCF = "DNA_VCF"
c_I_NUMBER_RETURNED_CLASS_ERRORS = 10
c_STR_NO_CALL = "NA"

# Constants for tab files
c_I_TAB_DNA_LOCATION = 0
c_I_TAB_DNA_REF = 1
c_I_TAB_DNA_CALL = 2
c_I_TAB_DNA_COVERAGE = 3
c_I_TAB_RNA_LOCATION = 4
c_I_TAB_RNA_REF = 5
c_I_TAB_RNA_CALL = 6
c_I_TAB_RNA_COVERAGE = 7
c_I_MIN_COVERAGE = 10

def func_make_weighted_coverage( str_RNA_coverage, str_DNA_coverage ):
  # Ad hoc waiting for coverage weight coverage by both amount of coverage and balance of coverage
  # ( prefering balanced maginitude of covarege between DNA and RNA )
  i_RNA_coverage = int( str_RNA_coverage ) * 1.0
  i_DNA_coverage = int( str_DNA_coverage ) * 1.0
  i_max = max( i_RNA_coverage, i_DNA_coverage )
  i_min = min( i_RNA_coverage, i_DNA_coverage )
  return ( 1 - ( ( i_max - i_min ) /  i_max ) ) * ( i_RNA_coverage + i_DNA_coverage )
#  return ( 1 - abs( max( min( 0, i_RNA_coverage / i_DNA_coverage ), 2 ) - 1 ) ) * ( i_RNA_coverage + i_DNA_coverage )

# Parse arguments
prsr_arguments = argparse.ArgumentParser( prog = "make_inspector_json.py", description = "Creates the json object needed to view a RNA-Seq mutation validation comparison run.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--input_files", required = True, dest = "lstr_input_files", action = "append", help = "A list of Sample name, RNA bam file, DNA bam file, RNA VCF file, DNA VCF file, tab file in a comma delimited string. Can be used more than once." )
prsr_arguments.add_argument( "--output_file", required = True, dest = "str_output_file", action = "store", help = "File to store the json object." )
args_call = prsr_arguments.parse_args()

# Object holding file info to become a json object for inspector visualization
# { sample : { RNA: file_path, DNA: file_path, FN: "Chr1-234 (145)": { "Chr": "1", "Loc": "234", "Cov": "145", "Ref": "A", "Alt": "T", "Strand": "+" }, FP: ..., TP: ...} }
dict_inspector = {}

# Go through each set of input files / info
for str_info in args_call.lstr_input_files:
  # Break up the input file string into file tokens
  str_sample_name, str_rna_bam, str_dna_bam, str_rna_vcf, str_dna_vcf, str_tab = str_info.split(",")

  # Check to make sure the files exist.
#  for str_file in [ str_rna_bam, str_dna_bam, str_rna_vcf, str_dna_vcf, str_tab ]:
#    if not os.path.exists( str_file ):
#      print( "Error. The input file " + str_file + " does not exist. Skipping sample " + str_sample_name + "." )
#      continue

  # Add sample and file names
  dict_sample = {}
  dict_sample[ c_STR_INSPECTOR_DNA_BAM ] = str_dna_bam
  dict_sample[ c_STR_INSPECTOR_RNA_BAM ] = str_rna_bam
  dict_sample[ c_STR_INSPECTOR_DNA_VCF ] = str_dna_vcf
  dict_sample[ c_STR_INSPECTOR_RNA_VCF ] = str_rna_vcf

  # Read in the VCF file
  llstr_tp_total = []
  llstr_tp = []
  llstr_fp_total = []
  llstr_fp = []
  llstr_fn_total = []
  llstr_fn = []
  # Open tab file
  with open( str_tab, "r" ) as  hndl_tab:
    csv_reader = csv.reader( hndl_tab, delimiter = "\t" )
    for lstr_tokens in csv_reader:
      # Skip the comments
      if lstr_tokens[0][0] == "#":
        continue

       # Skip NA coverage in either set
 #     # Change NA coverage to 0 coverage
      if lstr_tokens[ c_I_TAB_RNA_COVERAGE ].lower() == "na":
        continue
#        lstr_tokens[ c_I_TAB_RNA_COVERAGE ] = 0
      if lstr_tokens[ c_I_TAB_DNA_COVERAGE ].lower() == "na":
        continue
#        lstr_tokens[ c_I_TAB_DNA_COVERAGE ] = 0

      # Skip under a certain coverage in either DNA or RNA Coverage
      if int( lstr_tokens[ c_I_TAB_RNA_COVERAGE ] ) < ( c_I_MIN_COVERAGE + 1 ):
        continue
      if int( lstr_tokens[ c_I_TAB_DNA_COVERAGE ] ) < ( c_I_MIN_COVERAGE + 1 ):
        continue

      # Sort into error classes.
      STR_DNA_CALL = lstr_tokens[ c_I_TAB_DNA_CALL ]
      STR_RNA_CALL = lstr_tokens[ c_I_TAB_RNA_CALL ]

      # Append potential variants
      if STR_DNA_CALL == c_STR_NO_CALL:
        # FP
        if not STR_RNA_CALL == c_STR_NO_CALL:
          llstr_fp_total.append( lstr_tokens )
      # FN
      elif STR_RNA_CALL == c_STR_NO_CALL:
        llstr_fn_total.append( lstr_tokens )
      # TP
      else:
        llstr_tp_total.append( lstr_tokens )

  # Highest coverage no balance
  # Highest coverage with balance
  # Get the top number of FP by coverage
  llstr_fp_total.sort( key=lambda x: int( x[ c_I_TAB_RNA_COVERAGE ] ), reverse=True )
  llstr_fp = llstr_fp_total[0:c_I_NUMBER_RETURNED_CLASS_ERRORS]
  llstr_fp_total.sort( key=lambda x: func_make_weighted_coverage(x[ c_I_TAB_RNA_COVERAGE ], x[ c_I_TAB_DNA_COVERAGE ] ), reverse=True )
  llstr_fp.extend( llstr_fp_total[ 0:c_I_NUMBER_RETURNED_CLASS_ERRORS ] )
  # Get the top number of TP by coverage
  llstr_tp_total.sort( key=lambda x: int( x[ c_I_TAB_RNA_COVERAGE ] ), reverse=True )
  llstr_tp = llstr_tp_total[0:c_I_NUMBER_RETURNED_CLASS_ERRORS]
  llstr_tp_total.sort( key=lambda x: func_make_weighted_coverage(x[ c_I_TAB_RNA_COVERAGE ], x[ c_I_TAB_DNA_COVERAGE ] ), reverse=True )
  llstr_tp.extend( llstr_tp_total[ 0:c_I_NUMBER_RETURNED_CLASS_ERRORS ] )
  # Get the top number of FN by coverage
  llstr_fn_total.sort( key=lambda x: int( x[ c_I_TAB_DNA_COVERAGE ] ), reverse=True )
  llstr_fn = llstr_fn_total[0:c_I_NUMBER_RETURNED_CLASS_ERRORS]
  llstr_fn_total.sort( key=lambda x: func_make_weighted_coverage(x[ c_I_TAB_RNA_COVERAGE ], x[ c_I_TAB_DNA_COVERAGE ] ), reverse=True )
  llstr_fn.extend( llstr_fn_total[0:c_I_NUMBER_RETURNED_CLASS_ERRORS] )

  # Add error class info
  # FP
  dict_fp = {}
  for lstr_fp in llstr_fp:
    lstr_alt = list( set( lstr_fp[ c_I_TAB_RNA_CALL ].split("/") ) )
    lstr_alt = [ str_base for str_base in lstr_alt if str_base not in [ lstr_fp[ c_I_TAB_RNA_REF ] ]]
    lstr_temp_chr_loc = lstr_fp[ c_I_TAB_RNA_LOCATION ].split( "--" )
    str_temp_chr = lstr_temp_chr_loc[ 0 ][ 3: ] if ( len( lstr_temp_chr_loc[ 0 ] ) > 3 ) and ( lstr_temp_chr_loc[0][ 0:3 ].lower() == "chr" ) else lstr_temp_chr_loc[ 0 ]
    str_temp_loc = lstr_temp_chr_loc[ 1 ]
    for str_alt_base in lstr_alt:
      dict_temp = { "Chr": str_temp_chr, "Loc": str_temp_loc, "Cov_dna": lstr_fp[ c_I_TAB_DNA_COVERAGE ],
                    "Cov": lstr_fp[ c_I_TAB_RNA_COVERAGE ], "Ref": lstr_fp[ c_I_TAB_RNA_REF ],
                    "Alt": str_alt_base, "Strand": "+" }
      dict_fp[ "-".join( [ "Chr"+str_temp_chr, str_temp_loc ] ) + " (" + str( lstr_fp[ c_I_TAB_RNA_COVERAGE ] ) + ")" ] = dict_temp
  dict_sample[ c_STR_INSPECTOR_FP ] = dict_fp
  # TP
  dict_tp = {}
  for lstr_tp in llstr_tp:
    lstr_alt = list( set( lstr_tp[ c_I_TAB_RNA_CALL ].split("/") ) )
    lstr_alt = [ str_base for str_base in lstr_alt if str_base not in [ lstr_tp[ c_I_TAB_RNA_REF ] ]]
    lstr_temp_chr_loc = lstr_tp[ c_I_TAB_RNA_LOCATION ].split( "--" )
    str_temp_chr = lstr_temp_chr_loc[ 0 ][ 3: ] if ( len( lstr_temp_chr_loc[ 0 ] ) > 3 ) and ( lstr_temp_chr_loc[0][ 0:3 ].lower() == "chr" ) else lstr_temp_chr_loc[ 0 ]
    str_temp_loc = lstr_temp_chr_loc[ 1 ]
    for str_alt_base in lstr_alt:
      dict_temp = { "Chr": str_temp_chr, "Loc": str_temp_loc, "Cov_dna": lstr_tp[ c_I_TAB_DNA_COVERAGE ],
                    "Cov": lstr_tp[ c_I_TAB_RNA_COVERAGE ], "Ref": lstr_tp[ c_I_TAB_RNA_REF ],
                    "Alt": str_alt_base, "Strand": "+" }
      dict_tp[ "-".join( [ "Chr"+str_temp_chr, str_temp_loc ] ) + " (" + str( lstr_tp[ c_I_TAB_RNA_COVERAGE ] ) + ")" ] = dict_temp
  dict_sample[ c_STR_INSPECTOR_TP ] = dict_tp
  # FN
  dict_fn = {}
  for lstr_fn in llstr_fn:
    lstr_alt = list( set( lstr_fn[ c_I_TAB_DNA_CALL ].split("/") ) )
    lstr_alt = [ str_base for str_base in lstr_alt if str_base not in [ lstr_tp[ c_I_TAB_DNA_REF ] ]]
    lstr_temp_chr_loc = lstr_fn[ c_I_TAB_DNA_LOCATION ].split( "--" )
    str_temp_chr = lstr_temp_chr_loc[ 0 ][ 3: ] if ( len( lstr_temp_chr_loc[ 0 ] ) > 3 ) and ( lstr_temp_chr_loc[0][ 0:3 ].lower() == "chr" ) else lstr_temp_chr_loc[ 0 ]
    str_temp_loc = lstr_temp_chr_loc[ 1 ]
    for str_alt_base in lstr_alt:
      dict_temp = { "Chr": str_temp_chr, "Loc": str_temp_loc, "Cov_dna": lstr_fn[ c_I_TAB_DNA_COVERAGE ],
                    "Cov": lstr_fn[ c_I_TAB_RNA_COVERAGE ], "Ref": lstr_fn[ c_I_TAB_RNA_REF ],
                    "Alt": str_alt_base, "Strand": "+" }
      dict_fn[ "-".join( [ "Chr"+str_temp_chr, str_temp_loc ] ) + " (" + str( lstr_fn[ c_I_TAB_DNA_COVERAGE ] ) + ")" ] = dict_temp
  dict_sample[ c_STR_INSPECTOR_FN ] = dict_fn
  dict_inspector[ str_sample_name ] = dict_sample

# Open handle and write json object to file
with open( args_call.str_output_file, "w" ) as hndl_output:
  hndl_output.write( json.dumps( dict_inspector, sort_keys=True, indent=2 ) )
