#!/usr/bin/env python

# VCF contants
c_I_CHR_INDEX=0
c_I_POS_INDEX=1
c_CHR_VCF_DELIM="\t"
c_CHR_COMMENT="#"

# libraries
import argparse
import csv
from collections import deque
import os

# Commandline parameters
prsr_arguments = argparse.ArgumentParser( prog = "filter_variant_clusters.py", description = "Filter clusters of mutations on vcf files.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--window", dest = "i_base_window", default=35, type=int, action = "store", help = "Number of consequetive bases used as a window when identifying clusters." )
prsr_arguments.add_argument( "--cluster", dest = "i_max_allowed_cluster", default=2, type=int, action = "store", help = "The Maximum number of mutations allowed in a cluster in a windows length." )
prsr_arguments.add_argument( dest = "str_vcf_file_in", action = "store", help = "Input vcf file (coordinate sorted)." )
prsr_arguments.add_argument( dest = "str_output_file", action = "store", help = "Output tab file of snp calls." )
args = prsr_arguments.parse_args()

class Window:

  def __init__( self, hndl_out, i_window_length, i_cluster_length ):
    """
    Initialize Window with a size, allowable size of clusters, output file to write features that pass.

    * hndl_out : Output file to write features to that pass the filter
               : File handle
    * i_window_length : The length of the moving window for filtering.
                      : String
    * i_cluster_length : The max number of features allowed in a window.
                       : String
    """
    # Caches of each VCF line ( feature ) and if the line should be deleted
    self.llstr_line_cache = deque( [] )
    self.lf_line_keep = deque( [] )

    # Window length
    self.i_window_length = i_window_length
    # Max cluster membership
    self.i_cluster_length = i_cluster_length

    # Output file
    self.hndl_output_file = hndl_out

    # Current state
    self.str_prev_chr = "0"
    self.i_max_window_position = 0

  def func_flush( self ):
    """
    Flush variants currently in window, potentially writing features if they pass the filter.
    """
    # Remove each line (writing them potentially to file)
    for i_line in xrange( len( self.llstr_line_cache ) ):
      self.func_remove_item()
    self.i_min_window_position = 0

  def func_reduce_to_window( self ):
    """
    Check the span of the underlying data and make sure they are in a window of the allowable size.
    Remove any feature out which has fallen out of the window and potentially write to file.
    """
    # Update min position
    self.func_update_max_pos()
    i_number_lines_to_delete = 0
    for lstr_line in self.llstr_line_cache:
      if ( ( self.i_max_window_position - int( lstr_line[ c_I_POS_INDEX ] ) ) + 1 ) > self.i_window_length:
        i_number_lines_to_delete = i_number_lines_to_delete + 1
      else:
        break
    for i_remove in xrange( i_number_lines_to_delete ):
      self.func_remove_item()

  def func_update_max_pos( self ):
    """
    Update the minimum position of the window given the underlying features.
    """
    if ( self.llstr_line_cache ) == 0:
      self.i_max_window_position = 0
    else:
      self.i_max_window_position = int( self.llstr_line_cache[ -1 ][ c_I_POS_INDEX ] )

  def func_add_line( self, lstr_new_line ):
    """
    Add a line ot the window. Potentially write features that fall out of the span of the window.
    * lstr_new_line : List of strings representing a VCF feature.
                    : List of strings
    """

    # Current locations
    str_cur_chr = lstr_line[ c_I_CHR_INDEX ]

    # If new chromosome
    # Reset by throwing out all current mutations
    if not str_cur_chr == self.str_prev_chr:

      # Update chromosome
      self.str_prev_chr = str_cur_chr

      # Remove each line (writing them potentially to file)
      self.func_flush()

    # Add new item
    self.llstr_line_cache.append( lstr_line )
    self.lf_line_keep.append( True )

    # Remove items not within the window given the last item entered
    self.func_reduce_to_window()

    # If there are too many features, indicate that none of them can be kept
    if len( self.llstr_line_cache ) > self.i_cluster_length:
      self.lf_line_keep = deque( [ False for f_keep in self.lf_line_keep ] )

  def func_remove_item( self ):
    """
    Remove features from the window, writing them to a file if they pass the filter.
    """
    lstr_cur_line = self.llstr_line_cache.popleft()
    f_cur_keep = self.lf_line_keep.popleft()
    if f_cur_keep:
      self.hndl_output_file.write( c_CHR_VCF_DELIM.join( lstr_cur_line ) + "\n" )

# Window controlling writting to the output file
with open( args.str_output_file, "w" ) as hndl_out:
  window_current = Window( hndl_out, args.i_base_window, args.i_max_allowed_cluster )

  # Open input VCF file
  with open( args.str_vcf_file_in, "r" ) as hndl_vcf_in:
    # Read in line
    for lstr_line in csv.reader( hndl_vcf_in, delimiter = c_CHR_VCF_DELIM ):
      # Write comments directly to output
      if lstr_line[0][0] == c_CHR_COMMENT:
        hndl_out.write( c_CHR_VCF_DELIM.join( lstr_line ) + "\n" )
        continue

      # Add a feature into the window
      window_current.func_add_line( lstr_line )

  # Flush out / potentially write any remaining line in the window.
  window_current.func_flush()
