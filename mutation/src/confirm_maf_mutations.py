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
c_I_VCF_REF = 3
c_I_VCF_ALT = 4
c_I_VCF_GENOTYPE = 9

# MAF
c_I_MAF_CHR = 4
c_I_MAF_GENE_NAME = 0
c_I_MAF_START = 5
c_I_MAF_SAMPLE = 15
c_I_MAF_REF = 10
c_I_MAF_ALLELE_1 = 11
c_I_MAF_ALLELE_2 = 12

######## functions
def func_read_maf( str_maf_file ):
  
  # Dict
  # { sample: [ [ pos, gene ], [ pos, gene ] , ... ] }
  dict_maf = {}

  # Skip first line header
  f_skip = True

  with open( str_maf_file, "r" )  as hndl_maf:
    for str_line in csv.reader( hndl_maf, delimiter = "\t" ):
      if f_skip:
        f_skip = False
        continue
      str_pos = "--".join( [ "chr"+str_line[ c_I_MAF_CHR ], str_line[ c_I_MAF_START ] ] )
      str_sample = str_line[ c_I_MAF_SAMPLE ]
      str_gene_name = str_line[ c_I_MAF_GENE_NAME ].lower()
      str_ref = str_line[ c_I_MAF_REF ]
      str_genotype = "/".join( [ str_line[ c_I_MAF_ALLELE_1 ], str_line[ c_I_MAF_ALLELE_2 ] ] )
      dict_maf.setdefault( str_sample, [] ).append( [ str_pos, str_gene_name, str_ref, str_genotype ] )
  return dict_maf


def func_read_sample_mappings( str_file_sample_mappings ):
  
  dict_sample_mappings = {}
  f_skip = True

  with open( str_file_sample_mappings, "r" ) as hndl_mappings:
    for lstr_line in csv.reader( hndl_mappings, delimiter="\t" ):
      if f_skip:
        f_skip = False
        continue
      dict_sample_mappings[ lstr_line[ 0 ] ] = lstr_line[ 1 ]   
  return dict_sample_mappings


def func_read_vcf( str_file ):
  """
  Parse a VCF file.
  """

  # Dict to hold genomic locations
  # { chr--pos: None }
  dict_location = {}

  # Read in vcf file
  with open( str_file, "r" ) as hndl_vcf:
    for lstr_vcf_line in csv.reader( hndl_vcf, delimiter = "\t" ):
      # Skip comments
      if lstr_vcf_line[ 0 ][ 0 ] == "#":
        continue

      # Get genomic location
      lstr_ref_alt = [ lstr_vcf_line[ c_I_VCF_REF ] ] + lstr_vcf_line[ c_I_VCF_ALT ].split( "," )
      lstr_genotype = [ lstr_ref_alt[ int( str_index ) ] for str_index in lstr_vcf_line[ c_I_VCF_GENOTYPE ].split( ":" )[ 0 ].split( "/" ) ]
      str_genotype = "/".join( lstr_genotype )
      dict_location[ "--".join( [ lstr_vcf_line[ c_I_VCF_CHR ], lstr_vcf_line[ c_I_VCF_POS] ] ) ] = [ lstr_ref_alt[ 0 ], str_genotype ]
  return dict_location
      

prsr_arguments = argparse.ArgumentParser( prog = "confirm_maf_mutations.py", description = "Confirms mutations in maf file.", conflict_handler="resolve", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-m", "--maf", metavar = "Maf_file", dest = "str_maf", default = None, required = True, help = "Maf file." )
prsr_arguments.add_argument( "-s", "--sample", metavar = "Sample_mapping_file", dest = "str_sample_mappings", default = None, required = True, help = "Sample Mapping file." )
prsr_arguments.add_argument( "-k", "--key_genes", metavar = "Key_mutations", dest = "str_key_mappings", default = None, help = "Comma delimited gene name list. If given analysis will be reduced to these genes only.." )
prsr_arguments.add_argument( dest = "str_output_file", default = None, help = "PDF figure and related text file." )
args_call = prsr_arguments.parse_args()

print( "output" )
print( args_call.str_output_file )

# Holds the locations that were confirmed
# { gene : i_confirmed_instances }
dict_gene_confirmed = {}
# Holds the locations that were missed
# { gene : i_missed_instances }
dict_gene_missed = {}
# Holds the locations that had the wrong ref or genotype
# { str_pos : i_missed_ref/genotype  }
dict_wrong_genotype_pos = {}
# Holds the missed postions
llstr_missed_detail = []
# Samples not matched
lstr_samples_not_matched = []
# Number of vcf smaples that match
i_vcf_matched = 0

# Key genes (None is not being used)
lstr_key_genes = None
if args_call.str_key_mappings:
  lstr_key_genes = [ str_gene.lower() for str_gene in args_call.str_key_mappings.split(",") ]
  lstr_key_genes.reverse()

###### Make sure the same genotype.

# Read in maf file
# { sample : [ pos, gene ], [ pos, gene ], ... }
dict_maf_data = func_read_maf( args_call.str_maf )

# Read in sample mappings
#{ maf_sample : vcf_sample }
dict_sample_mappings = func_read_sample_mappings( args_call.str_sample_mappings )

# For vcf in maf file
for str_sample, llstr_positions in dict_maf_data.items():
  if str_sample in dict_sample_mappings:
    
    ## Get VCF SNPs
    print "Reading " + dict_sample_mappings[ str_sample ]
    dict_vcf_snps = func_read_vcf( dict_sample_mappings[ str_sample ] )
    i_vcf_matched = i_vcf_matched + 1
 
    # Check location
    for lstr_maf_pos in llstr_positions:
      # Make sure the gene is in the key genes
      if lstr_key_genes:
        if not lstr_maf_pos[ 1 ] in lstr_key_genes:
          continue

      if lstr_maf_pos[ 0 ] in dict_vcf_snps:
        lstr_cur_vcf_entry = dict_vcf_snps[ lstr_maf_pos[ 0 ] ]
        # Check correct genotype
        if lstr_maf_pos[ 2 ] == lstr_cur_vcf_entry[ 0 ] and lstr_maf_pos[ 3 ] == lstr_cur_vcf_entry[ 1 ]:
          dict_gene_confirmed[ lstr_maf_pos[1] ] = dict_gene_confirmed.setdefault( lstr_maf_pos[1], 0 ) + 1
        else:
          # { Position : [ [ mafref maf/genotype ], [ vcfref vcf/genotype ] ] }
          dict_wrong_genotype_pos[ lstr_maf_pos[ 0 ] ] = [ " ".join( [ lstr_maf_pos[ 2 ],lstr_maf_pos[ 3 ] ] ),
                                                           " ".join( [ lstr_cur_vcf_entry[ 0 ],
                                                           lstr_cur_vcf_entry[ 1 ] ] ), 
                                                           str_sample, 
                                                           lstr_maf_pos[ 1 ] ]
      else:
        dict_gene_missed[ lstr_maf_pos[1] ] = dict_gene_missed.setdefault( lstr_maf_pos[1], 0 ) + 1
        llstr_missed_detail.append( [ str_sample, dict_sample_mappings[ str_sample ], lstr_maf_pos[0], lstr_maf_pos[1] ] )
  else:
    lstr_samples_not_matched.append( str_sample )

# Output info
print "Count of VCF samples matched. ( " + str( i_vcf_matched ) + " ) "
print "Samples not matched in maf. ( " + str( len( set( lstr_samples_not_matched ) ) ) + " samples )"
#print "\n".join( set( lstr_samples_not_matched ) )
for str_gene in set( dict_gene_confirmed.keys() + dict_gene_missed.keys() ):
  print "Gene: " + str_gene + " Confirmed: " + str( dict_gene_confirmed.setdefault( str_gene, 0 ) ) + " Missed: " + str( dict_gene_missed.setdefault( str_gene, 0 ) )
if len( dict_wrong_genotype_pos ) > 0:
  print "Positions that had the wrong genotype:"
  for str_pos in dict_wrong_genotype_pos:
    lstr_cur_genotypes = dict_wrong_genotype_pos[ str_pos ]
    print "Sample: " + str_sample + " Position: " + str_pos + " Maf: " + lstr_cur_genotypes[ 0 ] + " VCF: " + lstr_cur_genotypes[ 1 ]

# Plot as a combined bar plot
lstr_write_genes = set( dict_gene_confirmed.keys() + dict_gene_missed.keys() )
if lstr_key_genes:
  lstr_write_genes = lstr_key_genes 
li_genes_indicies = [ ( i_index * 3 ) + 1.5 for i_index in range( len( lstr_key_genes ) ) ]
li_counts = []
i_cur_count_index = 0
li_counts_indices = []
lstr_count_colors = []
for str_gene in lstr_write_genes:
  li_counts.append( dict_gene_confirmed.get( str_gene, 0 ) )
  li_counts.append( dict_gene_missed.get( str_gene, 0 ) )
  i_cur_count_index = i_cur_count_index + 1
  li_counts_indices.append( i_cur_count_index )
  i_cur_count_index = i_cur_count_index + 1
  li_counts_indices.append( i_cur_count_index )
  i_cur_count_index = i_cur_count_index + 1
  lstr_count_colors.extend( [ "cyan", "magenta" ] )

fig, ax = plt.subplots()
li_bar_indices = range( 1,len( lstr_write_genes ) + 3 )
ax.barh( li_counts_indices, li_counts, color=lstr_count_colors )
plt.yticks( li_genes_indicies, lstr_write_genes )
ax.set_title( "SNPs in Exome Seq Recovered in RNA-Seq" )
ax.set_xlabel( "Count" )
ax.set_ylabel( "Gene" )
for tick in ax.yaxis.get_major_ticks():
  tick.label.set_fontsize(6)
plt.savefig( args_call.str_output_file )
plt.close( fig )
