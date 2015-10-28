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
import matplotlib.pyplot as plt
import os

# VCF related constants
c_I_VCF_CHR = 0
c_I_VCF_POS = 1
c_I_VCF_FILTER = 6
c_I_VCF_ANNOTATIONS = 7

# FLAGS to turn on optional filters
f_COMMON_FILTER = False
f_COVERAGE_FILTER = True
i_INT_MIN_COVERAGE = 10
f_GENE_MUTATION_COUNT_FILTER = True
i_MIN_GENE_MUTATION_COUNT = 8

i_gtf_annotations = 8
i_gtf_location = 0
i_tab_location = 0
i_tab_mutation = 1
i_gtf_mutation = 1
i_gtf_start_index = 3
i_gtf_stop_index = 4

# Tab related
# Number of columns between the first and second evidence group in the TAB
# (between DNA and RNA columns of the same data)
I_TAB_SECOND_SHIFT = 4

######## functions
def func_get_key_gene_locations( str_gtf_file, lstr_target_genes ):
  """
      Reads a GTF file and returns a dict of locations associate with a list of key gene names.
      { chr--loc: str_gene_name }
  """

  print "func_get_key_gene_locations"
  # Read the GTF file
  # Get the locations from the GFT file that match the gene names in the key mutations list
  dict_target_genes = {}
  with open( str_gtf_file, 'r' ) as hndl_gtf:
    for lstr_line in csv.reader( hndl_gtf, delimiter = "\t" ):
      lstr_annotations = lstr_line[ i_gtf_annotations ].split(" ")
      i_index = lstr_annotations.index( "gene_name" ) + 1
      str_gene_name = lstr_annotations[ i_index ].lower()[:-1].strip("\"")
      str_chr = lstr_line[ i_gtf_location ]
      # If in target genes keep
      if str_gene_name in lstr_target_genes:
        for int_position in range( int( lstr_line[ i_gtf_start_index ] ), int( lstr_line[ i_gtf_stop_index ] ) + 1 ):
          dict_target_genes[ "--".join( [ str_chr, str( int_position ) ] ) ] = str_gene_name
  return dict_target_genes


def func_plot_one( str_file_output, lstr_write_genes, li_write_mutations, lstr_mutation_order ):
  """
  Plot a single file.
  """

  print "func_plot_one"
  # Get order of key mutations
  lstr_gene_names = [ str_name.split(" " )[0] for str_name in lstr_write_genes ]
  li_index = [ lstr_mutation_order.index( str_key_name ) for str_key_name in lstr_gene_names ]

  # Reorder names and values
  lstr_write_genes = [ lstr_write_genes[ li_index.index( i_index ) ] for i_index in range( len( lstr_mutation_order ) ) ]
  lstr_write_genes.reverse()
  li_write_mutations = [ li_write_mutations[ li_index.index( i_index ) ] for i_index in range( len( lstr_mutation_order ) ) ]
  li_write_mutations.reverse()

  # Plot figure
  fig, ax = plt.subplots()
  li_bar_indices = range( 1,len( lstr_write_genes ) + 3 )
  ax.barh( li_bar_indices, [0] + li_write_mutations + [0] )
  plt.yticks( li_bar_indices, [""] + lstr_write_genes + [""] )
  ax.set_title( "Prevalence of Mutations in Study" )
  ax.set_ylabel( "Gene" )
  ax.set_xlabel( "Occurence (Count in Samples)" )
  plt.savefig( str_file_output )
  plt.close( fig )


def func_read_tab( str_file, lstr_target_genes, dict_target_genes_tab ):
  """
  Read in a tab file.
  """

  print "func_read_tab"
  # Initialize mutation counts mutation counts
  # { str_gene_name: 0 }
  dict_mutation_counts = dict( [ [ str_key_gene_name, 0 ] for str_key_gene_name in lstr_target_genes ] )
  # { str_gene_name: [] }
  dict_mutation_samples = dict( [ [ str_key_gene_name, [] ] for str_key_gene_name in lstr_target_genes ] )

  # Read in tab file
  with open( str_file, "r" ) as hndl_tab:
    for lstr_tab_line in csv.reader( hndl_tab, delimiter = "\t" ):
      # Make sure the entry is not NA (only in the tab file because there is an entry in the comparison sample
      if lstr_tab_line[ i_tab_location ] in dict_target_genes_tab:
        if lstr_tab_line[ i_tab_mutation ].lower() in [ "na", "0" ]:
          pass
        # Get the associated gene for the feature
        str_cur_gene = dict_target_genes_tabs[ lstr_tab_line[ i_tab_location ] ]
        # Update the gene's mutation count
        dict_mutation_counts[ str_cur_gene ] = dict_mutation_counts[ str_cur_gene ] + 1
        # Update that the sample had the mutation in it
        lstr_samples_cur_gene = dict_mutation_samples[ str_cur_gene ]
        if not str_file in lstr_samples_cur_gene:
          lstr_samples_cur_gene.append( str_file )
    return( ( dict_mutation_samples, dict_mutation_counts ) )


def func_read_vcf( str_file, lstr_target_genes, dict_target_genes_vcf ):
  """
  Parse a VCF file.
  """

  print "func_read_vcf"
  # Initialize mutation counts mutation counts
  # { str_gene_name: 0 }
  dict_mutation_counts = dict( [ [ str_key_gene_name, 0 ] for str_key_gene_name in lstr_target_genes ] )
  # { str_gene_name: [] }
  dict_mutation_samples = dict( [ [ str_key_gene_name, [] ] for str_key_gene_name in lstr_target_genes ] )

  # Read in vcf file
  with open( str_file, "r" ) as hndl_vcf:
    for lstr_vcf_line in csv.reader( hndl_vcf, delimiter = "\t" ):
      # Skip comments
      if lstr_vcf_line[ 0 ][ 0 ] == "#":
        continue

      # Get genomic location
      str_position = "--".join( [ lstr_vcf_line[ c_I_VCF_CHR ], lstr_vcf_line[ c_I_VCF_POS] ] )
      # If the location is in a gene of interest
      if str_position in dict_target_genes_vcf:

        # Skip SNPs that failed a filter if a filter was ran
        if not lstr_vcf_line[ c_I_VCF_FILTER ].lower() in [ "pass", "." ]:
          continue

        # Optional Common filter.
        if f_COMMON_FILTER:
          if "COMMON=1" in lstr_vcf_line[ c_I_VCF_ANNOTATIONS ]:
            continue

        # Optional 
        if f_COVERAGE_FILTER:
          lstr_coverage_tokens = [ str_token for str_token in lstr_vcf_line[ c_I_VCF_ANNOTATIONS ].split(";") if "DP=" in str_token ]
          if len( lstr_coverage_tokens ) > 1:
            # Should not find mutliple entries
            print "Warning found more than one DP entry."
            continue
          elif len( lstr_coverage_tokens ) == 1:
            # Found one annotation for DP so filter
            i_DP = int( lstr_coverage_tokens[ 0 ].split("=")[ 1 ] )
            if i_DP < i_INT_MIN_COVERAGE:
              continue
          else:
            # Filter out; no evidence of coverage
            continue

        # Get the associated gene for the feature
        str_cur_gene = dict_target_genes_vcf[ str_position ]
        # Update the gene's mutation count
        dict_mutation_counts[ str_cur_gene ] = dict_mutation_counts[ str_cur_gene ] + 1
        # Update that the sample had the mutation in it
        lstr_samples_cur_gene = dict_mutation_samples[ str_cur_gene ]
        if not str_file in lstr_samples_cur_gene:
          lstr_samples_cur_gene.append( str_file )

    #Optional filter
    if f_GENE_MUTATION_COUNT_FILTER:
      dict_mutation_samples_filtered = {}
      dict_mutation_counts_filtered = {}
      for str_key,i_count in dict_mutation_counts.items():
        if i_count >= i_MIN_GENE_MUTATION_COUNT:
          dict_mutation_samples_filtered[ str_key ] = i_count
        else:
          if str_key in dict_mutation_samples:
            if str_key in dict_mutation_samples[ str_key ]:
              dict_mutation_samples[ str_key ].remove( str_key )  
        dict_mutation_samples_filtered[ str_key ] = dict_mutation_samples[ str_key ] 
      dict_mutation_samples = dict_mutation_samples_filtered
      dict_mutation_counts = dict_mutation_counts_filtered
    return( ( dict_mutation_samples, dict_mutation_counts ) )


def func_write_data_to_file( str_data_file, dict_mutation_samples, dict_mutation_counts ):
  """
  Write details to file.
  """

  print "func_write_data_to_file"
  # Write to file
  lstr_write_genes = []
  li_write_mutations = []
  with open( str_data_file, "w" ) as hndl_out_txt:
    # Write Mutation info
    hndl_out_txt.write( "\t".join( [ "Gene", "Mutations", "Samples", "Percent_samples\n" ] ) )
    for str_key in dict_mutation_counts.keys():
      i_cur_sample_count = len( dict_mutation_samples[ str_key ] )
      i_percent_count = 0 if i_cur_sample_count == 0 else round( float( i_cur_sample_count ) / i_samples, 2 )
      hndl_out_txt.write( "\t".join( [ str_key, str( dict_mutation_counts[ str_key ] ), str( i_cur_sample_count ), str( i_percent_count ) + "\n" ] ) )
      lstr_write_genes.append( str_key + " (" + str( i_percent_count ) + ")" )
      li_write_mutations.append( i_cur_sample_count )
    # Write sample info
    hndl_out_txt.write( "\n\n########## Sample detail\n" )
    if args_call.f_second_entry:
      hndl_out_txt.write( "### Second Type\n" )
    else:
      hndl_out_txt.write( "### First Type\n" )
    for str_key, lstr_values in dict_mutation_samples.items():
      hndl_out_txt.write( "\t".join( [ str_key ] + [ "\t".join( [ os.path.basename( str_cur_sample_name ) for str_cur_sample_name in lstr_values ] ) + "\n" ] ) )
    return ( lstr_write_genes, li_write_mutations )

prsr_arguments = argparse.ArgumentParser( prog = "tab_to_percent_mutations.py", description = "Counts mutation prevalence from tab files", conflict_handler="resolve", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-g", "--gtf", metavar = "GTF_file", dest = "str_gtf", default = None, required = True, help = "GTF file to get gene locations." )
prsr_arguments.add_argument( "-k", "--key", metavar = "Key_genes", dest = "str_key_genes", default = None, required = True, help = "Key genes of interest to count mutations in. comma delimited. Example: gene1,gene2,gene3,gene4." )
prsr_arguments.add_argument( "-t", "--file", metavar = "Input_file", dest = "lstr_tab_files", default = [], action = "append", help = "Tab or vcf files files to parse for mutations" )
prsr_arguments.add_argument( "-o", "--out_file", metavar = "Output_file", dest = "str_output_file", default = None, required = True, help = "PDF figure." )
prsr_arguments.add_argument( "-s", "--second", dest = "f_second_entry", default = False, action="store_true", help = "Switches between reading in the first or second evidence type in the tab files." )
args_call = prsr_arguments.parse_args()

print( "output" )
print( args_call.str_output_file )

# Require lstr_tab files or lstr_vcf_files
# These are the files that are evaluated
if len( args_call.lstr_tab_files ) == 0:
  print "Please provide either tab files or vcf files to parse."
  exit( 99 )

# Shift reading if looking at the second entry in the tab file.
if args_call.f_second_entry:
  i_tab_location = i_tab_location + I_TAB_SECOND_SHIFT
  i_tab_mutation = i_tab_mutation + I_TAB_SECOND_SHIFT

# Split the list of key mutations and normalize to lower case
lstr_target_genes = [ str_target.lower() for str_target in args_call.str_key_genes.split( "," ) ]

# Text file for output
str_output_text_file = os.path.splitext( args_call.str_output_file )[ 0 ] + ".txt"

# Number of samples
i_samples = len( args_call.lstr_tab_files )

# Read the GTF file
# Get the locations from the GFT file that match the gene names in the key mutations list
# dict[ chr--loc ] = gene_name
dict_target_genes = func_get_key_gene_locations( str_gtf_file = args_call.str_gtf, lstr_target_genes = lstr_target_genes )

# Read in tab file
# { str_gene_name : 0 }
dict_mutation_counts_total = dict( [ [ str_key_gene_name, 0 ] for str_key_gene_name in lstr_target_genes ] )
# { str_gene_name : [] }
dict_mutation_samples_total = dict( [ [ str_key_gene_name, [] ] for str_key_gene_name in lstr_target_genes ] )

# Count mutations from each file.
f_single_plots = len( args_call.lstr_tab_files ) < 2
for str_file in args_call.lstr_tab_files:
  # Print cur file
  print "Reading input file:" + str_file
  # Per file view output file name
  str_per_file_text_file = os.path.splitext( str_output_text_file )[ 0 ] + "_" + os.path.basename( str_file ) + ".txt"
  str_per_file_pdf_file = os.path.splitext( str_per_file_text_file )[ 0 ] + ".pdf"
  print str_per_file_text_file
  print str_per_file_pdf_file

  # Read table and count key mutations
  dict_mutation_samples, dict_mutation_counts = func_read_vcf( str_file=str_file, lstr_target_genes=lstr_target_genes, dict_target_genes_vcf=dict_target_genes ) if os.path.splitext( str_file )[ 1 ] == ".vcf" else func_read_tab( str_file=str_file, lstr_target_genes=lstr_target_genes, dict_target_genes_tab=dict_taget_genes )

  print dict_mutation_samples
  print dict_mutation_counts
  if f_single_plots:
    # Write to file
    lstr_write_genes, li_write_mutations = func_write_data_to_file( str_data_file = str_per_file_text_file, 
                                                                  dict_mutation_samples = dict_mutation_samples, 
                                                                  dict_mutation_counts = dict_mutation_counts )

    # Plot figure
    func_plot_one( str_file_output = str_per_file_pdf_file,
                     lstr_write_genes = lstr_write_genes,
                     li_write_mutations = li_write_mutations,
                     lstr_mutation_order = lstr_target_genes )

  # Combine the results from the file with the global
  for str_key, str_value in dict_mutation_counts.items():
    dict_mutation_counts_total[ str_key ] = dict_mutation_counts[ str_key ] + dict_mutation_counts_total[ str_key ]
  for str_key, str_value in dict_mutation_samples.items():
    dict_mutation_samples_total[ str_key ].extend( dict_mutation_samples[ str_key ] )

# Write the combined data
str_total_text_file = os.path.splitext( str_output_text_file )[ 0 ] + "_total.txt"
str_total_pdf_file = os.path.splitext( str_total_text_file )[ 0 ] + ".pdf"
lstr_write_genes_total, li_write_mutations_total = func_write_data_to_file( str_data_file = str_total_text_file,
                                                                            dict_mutation_samples = dict_mutation_samples_total,
                                                                            dict_mutation_counts = dict_mutation_counts_total )

print "*********"
print str_total_pdf_file
print lstr_write_genes_total
print li_write_mutations_total
print lstr_target_genes 
# Plot the combined results
func_plot_one( str_file_output = str_total_pdf_file, lstr_write_genes = lstr_write_genes_total, li_write_mutations = li_write_mutations_total, lstr_mutation_order = lstr_target_genes )
