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

prsr_arguments = argparse.ArgumentParser( prog = "vcfs_to_genotype_matrix", description = "Collapses vcfs to genotype matrices of snps", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-l", "--list", dest = "str_output_list_file", default = "", action = "store", help = "Absolute path of output list file." )
prsr_arguments.add_argument( "-m", "--matrix", dest = "str_output_matrix_file", default = "", action = "store", help = "Absolute path of output matrix file." )
prsr_arguments.add_argument( "dirs", nargs = '+', help = "One or more directories of vcf files." )
args = prsr_arguments.parse_args()

# Holds samples info { str_sample : {location:genotype} }
dict_samples = dict()
set_rows = set()
lstr_columns = []

set_loc_genotype = set()

set_trans = set()

# Go through each directory
for str_dir in args.dirs:

    # Get all vcf files in the directory
    for str_file in glob.glob( str_dir+os.sep+"*.vcf" ):
        
        # Add row
        lstr_columns = lstr_columns + [ str_file ]
        
        # Reduce vcf files to matrix
        with open( str_file, 'r') as hndl_vcf:
            
                dict_samples[ str_file ] = dict()
                reader = csv.reader( hndl_vcf, delimiter = STR_VCF_DELIMITER )
                lstr_locs = []
                lstr_geno = []
                lstr_cur_loc_genotype = []

                for lstr_line in reader:
                
                    # Skip comments
                    if not lstr_line[0][0] == CHR_COMMENT:
                
                      # Skip if they do not pass
                      if lstr_line[ I_FILTER_INDEX ] == STR_PASS:
                          continue

                      # Get genotype
                      li_genotype = [ int( x ) for x in lstr_line[ I_GENOTYPE_INDEX ].split( CHR_ANNOT_DELIM )[ 0 ].split( CHR_GENOTYPE_DELIMITER ) ]
                      lstr_genotype = lstr_line[ I_REF_INDEX ].split( CHR_GENOTYPE_MUT ) + lstr_line[ I_ALT_INDEX ].split( CHR_GENOTYPE_MUT )
                
                      lstr_genotype = [ lstr_genotype[ i_index ] for i_index in li_genotype ]
                      # Check to make sure the genotype is SNP
                      if len( [ 1 for str_genotype in lstr_genotype if len( str_genotype ) > 1] ) > 0:
                          continue

                      str_location = lstr_line[ I_CHR_INDEX ] + "--" + lstr_line[ I_POS_INDEX ]
                      lstr_locs.append( str_location )
                      str_genotype = "->".join( lstr_genotype )
                      dict_samples[ str_file ][ str_location ] = str_genotype
                      lstr_geno.append( "->".join( lstr_genotype ) )
                      lstr_cur_loc_genotype.append( str_location + ".." + str_genotype + "\n" )

                set_rows = set_rows.union( set( lstr_locs ) )
                set_trans = set_trans.union( set( lstr_geno ) )
                set_loc_genotype = set_loc_genotype.union( set( lstr_cur_loc_genotype ) )
                lstr_cur_loc_genotype = []
set_rows = sorted( list( set_rows ) )

# Write the list file
if args.str_output_list_file:
    with open( args.str_output_list_file, "w" ) as hndl_vcf_list:
        for str_line in set_loc_genotype:
            hndl_vcf_list.write( str_line )

# Write to output file
if args.str_output_matrix_file:
    i_samples_count = len( lstr_columns )
    with open( args.str_output_matrix_file, "w" ) as hndl_out:
        writer = csv.writer( hndl_out, delimiter = STR_VCF_DELIMITER )
    
        # Write header
        writer.writerow( [ "" ] + lstr_columns )
    
        # For each location
        for str_loc in set_rows:
          # Write header
          writer.writerow( [ str_loc ] + [ dict_samples[ str_samples ].get( str_loc, "NA" ) for str_samples in lstr_columns ] )
