#!/usr/bin/env python

# Constants
CHR_COMMENT = "#"
STR_TAB_DELIMITER = "\t"
STR_CHROM = "Chromosome"
STR_CHROM_UPDATE = "CHROM"
STR_EMPTY_FILE = "No Data"
STR_POS = "Position"
STR_POS_UPDATE = "POS"
STR_CHASM_PVALUE = "CHASM p-value"
STR_CHASM_PVALUE_UPDATE = "CHASM_PVALUE"
STR_CHASM_FDR = "CHASM FDR"
STR_CHASM_FDR_UPDATE = "CHASM_FDR"
STR_VEST_PVALUE = "VEST p-value"
STR_VEST_PVALUE_UPDATE = "VEST_PVALUE"
STR_VEST_FDR = "VEST FDR"
STR_VEST_FDR_UPDATE = "VEST_FDR"

import argparse
import csv

prsr_arguments = argparse.ArgumentParser( prog = "groom_cravat_annotation.py", description = "Formats Cravat annotation files so they are easily used with bcftools. Also sorta by coordinate order.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "str_input_file", help = "Input tab file." )
prsr_arguments.add_argument( "str_output_file", help = "Output groomed tab file." )
args = prsr_arguments.parse_args()

# Read header indicator
f_header_read = False
i_chrom_index = -1
i_pos_index = -1
i_chasm_pvalue = -1
i_chasm_fdr = -1
i_vest_pvalue = -1
i_vest_fdr = -1
lstr_header_order = []

# Cache tab file so it can be coordinate sorted
lstr_comments = []
llstr_tab = []

def func_is_int( c_letter ):
  try:
    int( c_letter )
    return True
  except:
    return False

# Read in tab file
if args.str_input_file:
  with open( args.str_output_file, "w" ) as hndl_out:
    with open( args.str_input_file, "r" ) as hndl_in:
      for lstr_line in csv.reader( hndl_in, delimiter = STR_TAB_DELIMITER ):

        # Store blank lines (apart of the comment header)
        if not lstr_line:
          lstr_comments.append( STR_TAB_DELIMITER.join( lstr_line ))
          continue

        # Store comments
        if lstr_line[0][0] == CHR_COMMENT:
          lstr_comments.append( STR_TAB_DELIMITER.join( lstr_line ) )
          continue

        # Work with the body of the file
        # If the header has not been read yet, read.
        if not f_header_read:
          f_header_read = True
          i_chrom_index = lstr_line.index( STR_CHROM )
          i_pos_index = lstr_line.index( STR_POS )
          i_chasm_pvalue = lstr_line.index( STR_CHASM_PVALUE ) if STR_CHASM_PVALUE in lstr_line else -1
          i_chasm_fdr = lstr_line.index( STR_CHASM_FDR ) if STR_CHASM_FDR in lstr_line else -1
          i_vest_pvalue = lstr_line.index( STR_VEST_PVALUE ) if STR_VEST_PVALUE in lstr_line else -1
          i_vest_fdr = lstr_line.index( STR_VEST_FDR ) if STR_VEST_FDR in lstr_line else -1
          lstr_header_order = [ i_chrom_index, i_pos_index, i_chasm_pvalue, i_chasm_fdr, i_vest_pvalue, i_vest_fdr ]

          # Sort for sorting
          lstr_comments.append( STR_TAB_DELIMITER.join([ STR_CHROM_UPDATE, STR_POS_UPDATE, STR_CHASM_PVALUE_UPDATE,
                                                    STR_CHASM_FDR_UPDATE, STR_VEST_PVALUE_UPDATE, STR_VEST_FDR_UPDATE ]) )
          continue

        # Check to see if this is not a good run
        if lstr_line[ 0 ] == STR_EMPTY_FILE:
          llstr_tab.append( [ 0, 0, STR_TAB_DELIMITER.join( [ "NA" ] * len( lstr_header_order ) ) ] )
          break

        # Shuffle to header index and reduce
        # And store
        str_chrom = lstr_line[ i_chrom_index ]
        if not func_is_int( str_chrom[0] ):
          if str_chrom[0] in ["c","C"]:
            str_chrom = str_chrom[3:]
        if func_is_int( str_chrom[0] ):
          if len( str_chrom ) == 1:
            str_chrom = "0" + str_chrom
        i_pos = int( lstr_line[ i_pos_index ] )
        llstr_tab.append( [ str_chrom, i_pos, STR_TAB_DELIMITER.join( [ lstr_line[ i_index ] if i_index > 0 else "NA" for i_index in lstr_header_order ] ) ] )

    # Sort by chr and pos
    llstr_tab = sorted( llstr_tab, key = lambda x: ( x[0], x[1] ) )

    # Write to file
    hndl_out.write( "\n".join( lstr_comments ) )
    hndl_out.write( "\n" )
    hndl_out.write( "\n".join( [ lstr_tab_line[ 2 ] for lstr_tab_line in llstr_tab ] ) )
