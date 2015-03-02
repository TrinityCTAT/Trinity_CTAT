#!/usr/bin/env python

# Constant
STR_POS_SEP = "--"

# output file constants
I_NUM_EXAMPLES_TO_WRITE = 10
I_BASE_WINDOW = 20

# TAB constants
I_FIRST_ENTRY_DEPTH_INDEX = 3
I_SECOND_ENTRY_DEPTH_INDEX = 7
STR_SECOND_CHR_LOCATION_ENTRY = 4

import argparse
import csv
import os

prsr_arguments = argparse.ArgumentParser( prog = "tabs_comparison_to_igv_batch.py", description = "Selects several comparisons in tabbed file of comparisons and puts them in an IGV batch file to visually explore.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--maf_dna", required = True, dest = "str_maf_dna", action = "store", help = "Maf vs DNA." )
prsr_arguments.add_argument( "--maf_rna", required = True, dest = "str_maf_rna", action = "store", help = "Maf vs RNA." )
prsr_arguments.add_argument( "--dna_rna", required = True, dest = "str_dna_rna", action = "store", help = "DNA vs RNA." )
prsr_arguments.add_argument( "--maf_bam", required = True, dest = "str_maf_bam", action = "store", help = "Maf bam.")
prsr_arguments.add_argument( "--dna_bam", required = True, dest = "str_dna_bam", action = "store", help = "DNA bam.")
prsr_arguments.add_argument( "--rna_bam", required = True, dest = "str_rna_bam", action = "store", help = "RNA bam.")
prsr_arguments.add_argument( "--maf_dna_out", required = True, dest = "str_maf_dna_out", action = "store", help = "MAF DNA output file.")
prsr_arguments.add_argument( "--dna_not_rna_out", required = True, dest = "str_dna_not_rna_out", action = "store", help = "DNA not RNA output file.")
prsr_arguments.add_argument( "--maf_not_dna_out", required = True, dest = "str_maf_not_dna_out", action = "store", help = "MAF NOT DNA output file.")
prsr_arguments.add_argument( "--maf_rna_out", required = True, dest = "str_maf_rna_out", action = "store", help = "MAF RNA output file.")
prsr_arguments.add_argument( "--maf_not_rna_out", required = True, dest = "str_maf_not_rna_out", action = "store", help = "MAF NOT RNA output file.")
prsr_arguments.add_argument( "--dna_rna_not_maf_out", required = True, dest = "str_dna_rna_not_maf_out", action = "store", help = "DNA RNA not MAF output file.")
prsr_arguments.add_argument( "--maf_rna_not_dna_out", required = True, dest = "str_maf_rna_not_dna_out", action = "store", help = "MAF RNA NOT DNA output file.")
prsr_arguments.add_argument( "--maf_dna_rna_out", required = True, dest = "str_maf_dna_rna_out", action = "store", help = "MAF DNA RNA output file.")
prsr_arguments.add_argument( "--error_out", required = True, dest = "str_error_out", action = "store", help = "Error output file.")
prsr_arguments.add_argument( "--error_depth_out", required = True, dest = "str_error_depth_out", action = "store", help = "Error depth output file.")
args = prsr_arguments.parse_args()

# MAF and DNA
dict_maf_dna = {}
dict_maf_not_dna = {}
# MAF and RNA
dict_maf_rna = {}
dict_maf_not_rna = {}
# DNA and RNA
dict_dna_rna = {}
dict_dna_not_rna = {}

def func_batch_script( lstr_locations, str_output_directory, str_file_name, lstr_file_bams ):
  
  if len( lstr_locations ) < 1:
    return

  with open( os.path.join( str_output_directory, str_file_name ), "w" ) as hndl_out:
    # Load reference genome
    hndl_out.write( "new\ngenome hg19\n" )
    # Load bams
    for str_bam in lstr_file_bams:
      hndl_out.write( "load " + str_bam + "\n" )
    # Set output directory for pictures
    hndl_out.write( "snapshotDirectory " + str_output_directory + "\n" )
    # Take pictures
    for str_key in lstr_locations:
      # goto chr1:65,289,335-65,309,335
      str_chr, str_position = str_key.split( STR_POS_SEP )
      i_position = int( str_position )
      hndl_out.write( "goto " + str_chr + ":" + str( i_position - I_BASE_WINDOW ) + "-" + str( i_position + I_BASE_WINDOW ) + "\n" )
      hndl_out.write( "sort position\ncollapse\nsnapshot\n" )

def func_get_highest_read_depth( set_chr_pos, dict_read_depth ):
  ### Make list of lists to sort chr-pos by read depth
  llx_read_depth = [ [ str_pos, min( int( dict_read_depth[ str_pos ][ I_FIRST_ENTRY_DEPTH_INDEX ] ), int( dict_read_depth[ str_pos ][ I_SECOND_ENTRY_DEPTH_INDEX ] ) ) ] 
                       for str_pos in set_chr_pos if ( ( not dict_read_depth[ str_pos ][ I_FIRST_ENTRY_DEPTH_INDEX ] == "NA" ) and ( not dict_read_depth[ str_pos ][ I_SECOND_ENTRY_DEPTH_INDEX ] == "NA" ) ) ]
  ### Sort big to small (read depth)
  llx_read_depth.sort( key = lambda x: x[ 1 ], reverse=True )
  ### Reduce the list of lists to a list of chr-pos of the required length
  llx_read_depth = [ l_entry[0] for l_entry in llx_read_depth[ : min( I_NUM_EXAMPLES_TO_WRITE, len( llx_read_depth ) ) ] ]
  return( llx_read_depth )

# Read in the maf dna comparison file
with open( args.str_maf_dna, 'r' ) as hndle_maf_dna:
  print "Reading " + args.str_maf_dna
  reader_maf_dna = csv.reader( hndle_maf_dna, delimiter = "\t" )
  for lstr_maf_dna in reader_maf_dna:
    # TODO Quick fix to get scripts running check
    if lstr_maf_dna[ I_FIRST_ENTRY_DEPTH_INDEX -3 ] == "NA":
      continue
    if lstr_maf_dna[ STR_SECOND_CHR_LOCATION_ENTRY ] == "NA":
      dict_maf_not_dna[ lstr_maf_dna[ 0 ] ] = lstr_maf_dna
    else:
      dict_maf_dna[ lstr_maf_dna[ 0 ] ] = lstr_maf_dna

# Read in the maf rna comparison file
with open( args.str_maf_rna, 'r' ) as hndle_maf_rna:
  print "Reading " + args.str_maf_rna
  reader_maf_rna = csv.reader( hndle_maf_rna, delimiter = "\t" )
  for lstr_maf_rna in reader_maf_rna:
    # TODO Quick fix to get scripts runing check
    if lstr_maf_rna[ I_FIRST_ENTRY_DEPTH_INDEX - 3 ] == "NA":
      continue
    if lstr_maf_rna[ STR_SECOND_CHR_LOCATION_ENTRY ] == "NA":
      dict_maf_not_rna[ lstr_maf_rna[ 0 ] ] = lstr_maf_rna
    else:
      dict_maf_rna[ lstr_maf_rna[ 0 ] ] = lstr_maf_rna

# Read in the dna rna comparison file
with open( args.str_dna_rna, 'r' ) as hndle_dna_rna:
  print "Reading " + args.str_dna_rna
  reader_dna_rna = csv.reader( hndle_dna_rna, delimiter = "\t" )
  for lstr_dna_rna in reader_dna_rna:
    # TODO Quick fix to get scripts running check
    if lstr_dna_rna[ I_FIRST_ENTRY_DEPTH_INDEX - 3 ] == "NA":
      continue
    if lstr_dna_rna[ STR_SECOND_CHR_LOCATION_ENTRY ] == "NA":
      dict_dna_not_rna[ lstr_dna_rna[ 0 ] ] = lstr_dna_rna
    else:
      dict_dna_rna[ lstr_dna_rna[ 0 ] ] = lstr_dna_rna

# We want to look at the following
## In maf and in DNA
lstr_maf_dna_keys = dict_maf_dna.keys()
lstr_maf_dna_keys = lstr_maf_dna_keys[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_maf_dna_keys ) ) ]
func_batch_script( lstr_locations = lstr_maf_dna_keys, str_output_directory = os.path.split(args.str_maf_dna_out)[0], 
                   str_file_name = os.path.basename(args.str_maf_dna_out), lstr_file_bams = [ args.str_dna_bam ] )

## In maf and not in DNA
lstr_maf_not_dna_keys = dict_maf_not_dna.keys()
lstr_maf_not_dna_keys = lstr_maf_not_dna_keys[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_maf_not_dna_keys ) ) ]
func_batch_script( lstr_locations = lstr_maf_not_dna_keys, str_output_directory = os.path.split( args.str_maf_not_dna_out )[0], 
                   str_file_name = os.path.basename(args.str_maf_not_dna_out), lstr_file_bams = [ args.str_dna_bam ] )

## In maf and in RNA
lstr_maf_rna_keys = dict_maf_rna.keys()
lstr_maf_rna_keys = lstr_maf_rna_keys[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_maf_rna_keys ) ) ]
func_batch_script( lstr_locations = lstr_maf_rna_keys, str_output_directory = os.path.split( args.str_maf_rna_out )[0], 
                   str_file_name = os.path.basename(args.str_maf_rna_out), lstr_file_bams = [ args.str_dna_bam, args.str_rna_bam ] )

## In maf and not in RNA
lstr_maf_not_rna_keys = dict_maf_not_rna.keys()
lstr_maf_not_rna_keys = lstr_maf_not_rna_keys[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_maf_not_rna_keys ) ) ]
func_batch_script( lstr_locations = lstr_maf_not_rna_keys, str_output_directory = os.path.split(args.str_maf_not_rna_out)[0], 
                   str_file_name = os.path.basename(args.str_maf_not_rna_out), lstr_file_bams = [ args.str_dna_bam, args.str_rna_bam ] )

## In DNA and not in RNA
lstr_dna_not_rna_keys = dict_dna_not_rna.keys()
lstr_dna_not_rna_keys = lstr_dna_not_rna_keys[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_dna_not_rna_keys ) ) ]
func_batch_script( lstr_locations = lstr_dna_not_rna_keys, str_output_directory = os.path.split(args.str_dna_not_rna_out)[0], 
                   str_file_name = os.path.basename(args.str_dna_not_rna_out), lstr_file_bams = [ args.str_dna_bam, args.str_rna_bam ] )

## In DNA and RNA but not in MAF ( and high read depth )
set_str_dna_rna = set( dict_dna_rna.keys() )
set_str_maf_dna = set( dict_maf_dna.keys() )
### Get only in DNA_RNA
set_str_dna_rna_not_maf = set_str_dna_rna - set_str_maf_dna
lstr_pos = func_get_highest_read_depth( list( set_str_dna_rna_not_maf ), dict_dna_rna )
func_batch_script( lstr_locations = lstr_pos, str_output_directory = os.path.split( args.str_dna_rna_not_maf_out )[0], 
                   str_file_name = os.path.basename( args.str_dna_rna_not_maf_out ), lstr_file_bams = [ args.str_dna_bam, args.str_rna_bam ] )

## In MAF and RNA but not DNA ( and high read depth )
set_str_maf_rna = set( dict_maf_rna.keys() )
### Get things in MAF and RNA but not DNA
set_str_maf_rna_not_dna = set_str_maf_rna - set_str_dna_rna
lstr_pos = func_get_highest_read_depth( list( set_str_dna_rna ), dict_dna_rna )
func_batch_script( lstr_locations = lstr_pos, str_output_directory = os.path.split(args.str_maf_rna_not_dna_out)[0], 
                   str_file_name = os.path.basename(args.str_maf_rna_not_dna_out), lstr_file_bams = [ args.str_dna_bam, args.str_rna_bam ] )  

## In MAF and DNA and RNA
lstr_all_keys = list( set_str_dna_rna.intersection( set_str_maf_dna ) )
lstr_all_keys = lstr_all_keys[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_all_keys ) ) ]
func_batch_script( lstr_locations = lstr_all_keys, str_output_directory = os.path.split(args.str_maf_dna_rna_out)[0], 
                   str_file_name = os.path.basename(args.str_maf_dna_rna_out), lstr_file_bams = [ args.str_dna_bam, args.str_rna_bam ] )

## QC look for DNA seq entries that have calls but no read depth
lstr_dna_depth_error = [ str_dna_depth_error for str_dna_depth_error in set_str_dna_rna if dict_dna_rna[ str_dna_depth_error ][ I_FIRST_ENTRY_DEPTH_INDEX ] == "NA" ]
lstr_dna_depth_error = lstr_dna_depth_error[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_dna_depth_error ) ) ]
func_batch_script( lstr_locations = lstr_dna_depth_error, str_output_directory = os.path.split(args.str_error_out)[0], 
                   str_file_name = os.path.basename(args.str_error_out), lstr_file_bams = [ args.str_dna_bam ] )

## QC look for RNA seq entries that have calls but no read depth
lstr_rna_depth_error = [ str_rna_depth_error for str_rna_depth_error in set_str_dna_rna if ( ( not dict_dna_rna[ str_dna_depth_error ][ I_FIRST_ENTRY_DEPTH_INDEX + 1 ] == "NA" ) and ( dict_dna_rna[ str_dna_depth_error ][ I_SECOND_ENTRY_DEPTH_INDEX ] == "NA" ) ) ]
lstr_rna_depth_error = lstr_rna_depth_error[ : min( I_NUM_EXAMPLES_TO_WRITE, len( lstr_rna_depth_error ) ) ]
func_batch_script( lstr_locations = lstr_rna_depth_error, str_output_directory = os.path.split( args.str_error_depth_out)[0], 
                   str_file_name = os.path.basename(args.str_error_depth_out), lstr_file_bams = [ args.str_rna_bam ] )
