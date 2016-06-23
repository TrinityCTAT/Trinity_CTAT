#!/usr/bin/env python

# Constants
CHR_ANNOT_DELIM = ":"
CHR_COMMENT = "#"
CHR_GENOTYPE_DELIMITER = "/"
CHR_GENOTYPE_MUT = ","

STR_VCF_DELIMITER = "\t"
STR_PASS = "PASS"

I_CHR_INDEX = 0
I_POS_INDEX = 1
I_REF_INDEX = 3
I_ALT_INDEX = 4
I_FILTER_INDEX = 6
I_GENOTYPE_INDEX = 9

import argparse
import csv
import glob
import os

prsr_arguments = argparse.ArgumentParser( prog = "vcfs_to_genotype_matrix", description = "Collapses vcfs to genotype matrices", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-m", "--matrix", dest = "str_output_matrix_file", default = "", action = "store", help = "Absolute path of output matrix file." )
prsr_arguments.add_argument( "files", nargs = '+', help = "One or more directories of vcf files." )
prsr_arguments.add_argument( "--no_prefilter", action="store_true", dest="f_no_prefilter_mode", help="When the VCF file has not been filtered (indicated by using this flag), this script will not require the SNPs to pass a filter")
args = prsr_arguments.parse_args()

# Holds samples info { str_sample : {location:genotype} }
dict_samples = dict()
# The row features collected from the vcfs ( gentoypes )
set_rows = set()
# The columns ( samples )
lstr_columns = []

# Go through each directory
for str_file in args.files:
  # Add row
  lstr_columns = lstr_columns + [ os.path.basename( str_file ) ]

  # Reduce vcf files to matrix
  with open( str_file, 'r') as hndl_vcf:

    # Add a sample (current vf file)
    dict_samples[ os.path.basename( str_file ) ] = dict()
    # Create a file reader
    reader = csv.reader( hndl_vcf, delimiter = STR_VCF_DELIMITER )
    # Holds locations until they are added to rows
    lstr_locs = []

    # Read through VCF file
    for lstr_line in reader:

      # Skip comments
      if lstr_line[0][0] == CHR_COMMENT:
        continue

#        # Skip if they do not pass
#        if not args.f_no_prefilter_mode:
#          if not lstr_line[ I_FILTER_INDEX ] == STR_PASS:
#            continue

      # Get genotype
      li_genotype = [ int( x ) for x in lstr_line[ I_GENOTYPE_INDEX ].split( CHR_ANNOT_DELIM )[ 0 ].split( CHR_GENOTYPE_DELIMITER ) ]
      lstr_genotype = lstr_line[ I_REF_INDEX ].split( CHR_GENOTYPE_MUT ) + lstr_line[ I_ALT_INDEX ].split( CHR_GENOTYPE_MUT )
      lstr_genotype = [ lstr_genotype[ i_index ] for i_index in li_genotype ]

#        # Check to make sure the genotype is SNP or not a .
#        if len( [ 1 for str_genotype in lstr_genotype if ((len( str_genotype ) > 1) or (str_genotype == "."))] ) > 0:
#          continue

      # Make genomic location (chr--position) and store to row features
      str_location = lstr_line[ I_CHR_INDEX ] + "--" + lstr_line[ I_POS_INDEX ]
      lstr_locs.append( str_location )

      # Store genotype by location and file name
      str_genotype = "->".join( lstr_genotype )
      dict_samples[ os.path.basename( str_file ) ][ str_location ] = str_genotype

      # Update the matrix rows
      set_rows = set_rows.union( set( lstr_locs ) )

# Order the rows
set_rows = sorted( list( set_rows ) )

# Write to output file
if args.str_output_matrix_file:
  # get the number of columns (file for the matrix
  i_samples_count = len( lstr_columns )

  # Write file
  with open( args.str_output_matrix_file, "w" ) as hndl_out:

    writer = csv.writer( hndl_out, delimiter = STR_VCF_DELIMITER )

    # Write header
    writer.writerow( [ "" ] + lstr_columns )

    # For each location
    for str_loc in set_rows:
      # Write header
      writer.writerow( [ str_loc ] + [ dict_samples[ str_samples ].get( str_loc, "NA" ) for str_samples in lstr_columns ] )
