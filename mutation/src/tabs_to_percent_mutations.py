#!/usr/bin/env python

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2015"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

import argparse
import csv
import os

i_gtf_annotations = 8
i_gtf_location = 0
i_tab_location = 0
i_tab_mutation = 1
i_gtf_mutation = 1
i_gtf_start_index = 3
i_gtf_stop_index = 4

prsr_arguments = argparse.ArgumentParser( prog = "tab_to_percent_mutations.py", description = "Counts mutation prevalence from tab files", conflict_handler="resolve", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-g", "--gtf", metavar = "GTF_file", dest = "str_gtf", default = None, required = True, help = "GTF file to get gene locations." )
prsr_arguments.add_argument( "-k", "--key", metavar = "Key_genes", dest = "str_key_genes", default = None, required = True, help = "Key genes of interest to count mutations in. comma delimited. Example: gene1,gene2,gene3,gene4." )
prsr_arguments.add_argument( "-t", "--tab", metavar = "Tab_file", dest = "lstr_tab_files", default = None, required = True, action = "append", help = "Tab files to parse for mutations" )
prsr_arguments.add_argument( "-o", "--out_file", metavar = "Output_file", dest = "str_output_file", default = None, required = True, help = "PDF figure." )
prsr_arguments.add_argument( "-s", "--second", dest = "f_second_entry", default = False, action="store_true", help = "Switches between reading in the first or second evidence type in the tab files." )
args_call = prsr_arguments.parse_args()

# Shift reading if looking at the second entry in the tab file.
i_tab_second_shift = 4
if args_call.f_second_entry:
  i_tab_location = i_tab_location + i_tab_second_shift
  i_tab_mutation = i_tab_mutation + i_tab_second_shift

lstr_target_genes = [ str_target.lower() for str_target in args_call.str_key_genes.split( "," ) ]
# dict[ number ] = gene_name
dict_target_genes = {}
# Text file for output
str_output_text_file = os.path.splitext( args_call.str_output_file )[ 0 ] + ".txt"
# Number of samples
i_samples = 0

# Read the GTF file
with open( args_call.str_gtf, 'r' ) as hndl_gtf:
  for lstr_line in csv.reader( hndl_gtf, delimiter = "\t" ):
    lstr_annotations = lstr_line[ i_gtf_annotations ].split(" ")
    i_index = lstr_annotations.index( "gene_name" ) + 1
    str_gene_name = lstr_annotations[ i_index ].lower()[:-1].strip("\"")
    str_chr = lstr_line[ i_gtf_location ]

    if str_gene_name in lstr_target_genes:
      for int_position in range( int( lstr_line[ i_gtf_start_index ] ), int( lstr_line[ i_gtf_stop_index ] ) + 1 ):
        dict_target_genes[ "--".join( [ str_chr, str( int_position ) ] ) ] = str_gene_name

# Read in tab file
dict_mutation_counts = dict( [ [ str_key_gene_name, 0 ] for str_key_gene_name in lstr_target_genes ] )
dict_mutation_samples = dict( [ [ str_key_gene_name, [] ] for str_key_gene_name in lstr_target_genes ] )
for str_tab_file in args_call.lstr_tab_files:
  with open( str_tab_file, "r" ) as hndl_tab:
    i_samples = i_samples + 1
    for lstr_tab_line in csv.reader( hndl_tab, delimiter = "\t" ):
      if lstr_tab_line[ i_tab_location ] in dict_target_genes:
        if lstr_tab_line[ i_tab_mutation ].lower() == "na":
          pass
        str_cur_gene = dict_target_genes[ lstr_tab_line[ i_tab_location ] ]
        dict_mutation_counts[ str_cur_gene ] = dict_mutation_counts[ str_cur_gene ] + 1
        lstr_samples_cur_gene = dict_mutation_samples[ str_cur_gene ]
        if not str_tab_file in lstr_samples_cur_gene:
          lstr_samples_cur_gene.append( str_tab_file )

# Write to file
lstr_write_genes = []
li_write_mutations = []
with open( str_output_text_file, "w" ) as hndl_out_txt:
  # Write Mutation info
  hndl_out_txt.write( "\t".join( [ "Gene", "Mutations", "Samples", "Percent_samples\n" ] ) )
  for str_key in dict_mutation_counts.keys():
    i_cur_sample_count = len( dict_mutation_samples[ str_key ] )
    hndl_out_txt.write( "\t".join( [ str_key, str( dict_mutation_counts[ str_key ] ), str( i_cur_sample_count ), str( round( float( i_cur_sample_count ) / i_samples, 2 )) + "\n" ] ) )
    lstr_write_genes.append( str_key + " (" + str( round( i_cur_sample_count / i_samples, 2 ) ) + ")" )
    li_write_mutations.append( i_cur_sample_count )

  # Write sample info
  hndl_out_txt.write( "\n\n########## Sample detail\n" )
  if args_call.f_second_entry:
    hndl_out_txt.write( "### Second Type\n" )
  else:
    hndl_out_txt.write( "### First Type\n" )
  for str_key, lstr_values in dict_mutation_samples.items():
    hndl_out_txt.write( "\t".join( [ str_key ] + [ "\t".join( [ os.path.basename( str_cur_sample_name ) for str_cur_sample_name in lstr_values ] ) + "\n" ] ) )

# Plot figure
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
li_bar_indices = range( 1,len( lstr_write_genes ) + 3 )
ax.barh( li_bar_indices, [0] + li_write_mutations + [0] )
plt.yticks( li_bar_indices, [""] + lstr_write_genes + [""] )
ax.set_title( "Prevalence of Mutations in Study" )
ax.set_ylabel( "Gene" )
ax.set_xlabel( "Occurence (Count in Samples)" )
plt.savefig( args_call.str_output_file )
