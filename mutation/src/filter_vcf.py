#!/usr/bin/env python

# Constant
c_STR_POS_SEP = "--"
c_CHR_GENOTYPE_SEP = "/"

# VCF contants
c_I_CHR_INDEX = 0
c_I_POS_INDEX = 1
c_I_FILTER_INDEX = 6
c_CHR_VCF_DELIM = "\t"
c_STR_PASS = "PASS"
c_CHR_COMMENT = "#"
c_CHR_MONOMORPHIC_REFERENCE = "."

import argparse
import csv
import os

prsr_arguments = argparse.ArgumentParser( prog = "filter_vcf", description = "Perform filtering on vcf files.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-i", "--in", required = True, dest = "str_vcf_file_in", action = "store", help = "Input vcf file." )
prsr_arguments.add_argument( dest = "str_output_file", action = "store", help = "Output tab file of snp calls." )
args = prsr_arguments.parse_args()

# Current chr
str_current_chr = "chr0"
i_current_min_pos = 0

# Cache
llstr_line_cache = []
lf_line_keep = []


def func_add_line( llstr_lines, lstr_new_line, lf_keep_status ):

  llstr_line.append( lstr_line )
  lf_keep_status.append( True )
  # When adding check to make sure the number of features in the window have not gotten too big
  # If there are too many features, indicate that none of them can be kept
  if len( llstr_line ) > c_I_MAX_CLUSTER:
    lf_keep_status = [ False for f_keep in lf_keep_status ]

def func_remove_item( llstr_line, lf_keep_status, hndl_vcf_out ):

  lstr_cur_line = llstr_line.pop()
  f_cur_keep = lf_keep_status.pop()
  if f_cur_keep:
    hndl_vcf_out.write( c_CHR_VCF_DELIM.join( lstr_cur_line + ["\n"] ) )
  return integer( llstr_line[ 0 ][ c_I_POS_INDEX ] ) if len( llstr_line ) else 0      

# Open VCF file
print "Reading File: " + args.str_vcf_file_in
with open( args.str_vcf_file_in, "r" ) as hndl_vcf_in:
  with open( args.str_output_file, "w" ) as hndl_vcf_out:
    # Read in line
    for lstr_line in csv.reader( hndl_vcf_in, delimiter = c_CHR_VCF_DELIM ):
      # Write comments back
      if lstr_line[0][0] == c_CHR_COMMENT:
        hndl_vcf_out.write( c_CHR_VCF_DELIM.join( lstr_line + ["\n"] ) )

      # Current locations
      str_cur_chr = lstr_line[ [ c_I_CHR_INDEX ] ]
      i_cur_position = integer( lstr_line[ c_I_POS_INDEX ] )

      ## If new chromosome, remove all items
      if len( llstr_line_cache ) < c_I_MAX_CLUSTER:
        for i_line in xrange( len( llstr_line_cache ) ):
          func_remove_item( llstr_line=llstr_line_cache, lf_keep_status=lf_line_cache, hndl_vcf_out=hndl_vcf_out )
        i_current_min_pos = 0

      ## Add line in queue if it is within a window of c_I_WINDOW of the first entry
      if i_cur_position - i_current_min_pos <= c_I_WINDOW:
        # Add line
        func_add_line( llstr_lines=llstr_line_cache, lstr_new_line=lstr_line, lf_keep_status=lf_line_cache )
      else:
        # If not within the window
        ## Remove lines ( adding or not ) until the earlier entries are within the window.
        for lstr_line_new_win in llstr_line_cache:
          if i_cur_position - integer( lstr_line_new_win[ c_I_POS_INDEX ] ) > c_I_WINDOW:
            # If not within window, add line if appropriate and remove item
            i_cur_position = func_remove_item( llstr_line=llstr_line_cache, lf_keep_status=lf_line_cache, hndl_vcf_out=hndl_vcf_out )
          else:
            # If now within window, stop and add line
            func_add_line( llstr_lines=llstr_line_cache, lstr_new_line=lstr_line_new_win, lf_keep_status=lf_line_cache )
            break
