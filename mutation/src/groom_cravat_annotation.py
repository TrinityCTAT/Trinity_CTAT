#!/usr/bin/env python

# Constants
CHR_COMMENT = "#"
STR_TAB_DELIMITER = "\t"
STR_CHROM = "Chromosome"
STR_CHROM_UPDATE = "CHROM"
STR_POS = "Position"
STR_POS_UPDATE = "POS"
STR_CHASM_PVALUE = "CHASM cancer driver p-value (missense)"
STR_CHASM_PVALUE_UPDATE = "CHASM_PVALUE"
STR_CHASM_FDR = "CHASM cancer driver FDR (missense)"
STR_CHASM_FDR_UPDATE = "CHASM_FDR"
STR_VEST_PVALUE = "VEST pathogenicity p-value (non-silent)"
STR_VEST_PVALUE_UPDATE = "VEST_PVALUE"
STR_VEST_FDR = "VEST pathogenicity FDR (non-silent)"
STR_VEST_FDR_UPDATE = "VEST_FDR"

import argparse
import csv

prsr_arguments = argparse.ArgumentParser( prog = "groom_cravat_annotation.py", description = "Formats Cravat annotation files so they are easily used with bcftools.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "str_input_file", help = "Input tab file." )
prsr_arguments.add_argument( "str_output_file", help = "Output groomed tab file." )
args = prsr_arguments.parse_args()

# Stores the tab info
lstr_tab = []

# Read header indicator
f_header_read = False
i_chrom_index = -1
i_pos_index = -1
i_chasm_pvalue = -1
i_chasm_fdr = -1
i_vest_pvalue = -1
i_vest_fdr = -1
lstr_header_order = []

# Read in tab file
if args.str_input_file:
  with open( args.str_output_file, "w" ) as hndl_out:
    with open( args.str_input_file, "r" ) as hndl_in:
      for lstr_line in csv.reader( hndl_in, delimiter = STR_TAB_DELIMITER ):

        # Store blank lines
        if not lstr_line:
          lstr_tab.append( STR_TAB_DELIMITER.join( lstr_line ))
          continue

        # Store comments
        if lstr_line[0][0] == CHR_COMMENT:
          lstr_tab.append( STR_TAB_DELIMITER.join( lstr_line ) )
          continue

        # Work with the body of the file
        # If the header has not been read yet, read.
        if not f_header_read:
          f_header_read = True
          i_chrom_index = lstr_line.index( STR_CHROM )
          i_pos_index = lstr_line.index( STR_POS )
          i_chasm_pvalue = lstr_line.index( STR_CHASM_PVALUE )
          i_chasm_fdr = lstr_line.index( STR_CHASM_FDR )
          i_vest_pvalue = lstr_line.index( STR_VEST_PVALUE )
          i_vest_fdr = lstr_line.index( STR_VEST_FDR )
          lstr_header_order = [ i_chrom_index, i_pos_index, i_chasm_pvalue, i_chasm_fdr, i_vest_pvalue, i_vest_fdr ]
          lstr_tab.append( STR_TAB_DELIMITER.join([ STR_CHROM_UPDATE, STR_POS_UPDATE, STR_CHASM_PVALUE_UPDATE,
                                                    STR_CHASM_FDR_UPDATE, STR_VEST_PVALUE_UPDATE, STR_VEST_FDR_UPDATE ]) )

        # Shuffle to header index and reduce
        # And store
        lstr_tab.append( STR_TAB_DELIMITER.join( [ lstr_line[ i_index ] for i_index in lstr_header_order ] ) )

    for str_out_line in lstr_tab:
      hndl_out.write( str_out_line + "\n" )
