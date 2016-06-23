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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

c_I_MAF_RETURN_POS = 0
c_I_MAF_RETURN_GENE_NAME = 1
c_I_MAF_RETURN_REF = 2
c_I_MAF_RETURN_GENOTYPE = 3

c_I_VCF_RETURN_REF = 0
c_I_VCF_RETURN_GENOTYPE = 1

######## functions
def func_read_maf( str_maf_file ):
  """
  Read in a MAF file.
  Return a dict[ sample_name ] = [ [position info],[position info],...]
  where position info is [position, gene name, genotype]
  """
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
  """
  Read in a file that maps between MAF ids and samples (VCF files)
  This is a tab delimited file.
  MAF_Id\tVCF_file

  Return a dict[ MAF_id ] = VCF_name
  """
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

  Return a dict[ chr--pos ] = [ ALT, Genotype ]
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

# Holds the locations that were confirmed
# { gene : i_confirmed_instances }
dict_gene_confirmed = {}
# Holds the locations that were missed
# { gene : i_missed_instances }
dict_gene_missed = {}
# { gene : [ confirmed, sample,...] }
dict_gene_confirmed_sample = {}
# Holds the locations that were missed
# { gene : [ missed, sample,...] }
dict_gene_missed_sample = {}
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
        if not lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] in lstr_key_genes:
          continue

      if lstr_maf_pos[ c_I_MAF_RETURN_POS ] in dict_vcf_snps:
        lstr_cur_vcf_entry = dict_vcf_snps[ lstr_maf_pos[ c_I_MAF_RETURN_POS ] ]
        # Check correct genotype
        if lstr_maf_pos[ c_I_MAF_RETURN_REF ] == lstr_cur_vcf_entry[ c_I_VCF_RETURN_REF ] and lstr_maf_pos[ c_I_MAF_RETURN_GENOTYPE ] == lstr_cur_vcf_entry[ c_I_VCF_RETURN_GENOTYPE ]:
          dict_gene_confirmed[ lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] ] = dict_gene_confirmed.setdefault( lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ], 0 ) + 1
          dict_gene_confirmed_sample.setdefault( lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ], [] ).append( str_sample )
        else:
          # { Position : [ [ mafref maf/genotype ], [ vcfref vcf/genotype ] ] }
          dict_wrong_genotype_pos[ lstr_maf_pos[ c_I_MAF_RETURN_POS ] ] = [ " ".join( [ lstr_maf_pos[ c_I_MAF_RETURN_REF ],lstr_maf_pos[ c_I_MAF_RETURN_GENOTYPE ] ] ),
                                                           " ".join( [ lstr_cur_vcf_entry[ c_I_VCF_RETURN_REF ],
                                                           lstr_cur_vcf_entry[ c_I_VCF_RETURN_GENOTYPE ] ] ),
                                                           str_sample,
                                                           lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] ]
          dict_gene_missed[ lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] ] = dict_gene_missed.setdefault( lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ], 0 ) + 1
          llstr_missed_detail.append( [ str_sample, dict_sample_mappings[ str_sample ], lstr_maf_pos[ c_I_MAF_RETURN_POS ], lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] ] )
          dict_gene_missed_sample.setdefault( lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ], [] ).append( str_sample )
      else:
        dict_gene_missed[ lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] ] = dict_gene_missed.setdefault( lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ], 0 ) + 1
        llstr_missed_detail.append( [ str_sample, dict_sample_mappings[ str_sample ], lstr_maf_pos[ c_I_MAF_RETURN_POS ], lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ] ] )
        dict_gene_missed_sample.setdefault( lstr_maf_pos[ c_I_MAF_RETURN_GENE_NAME ], [] ).append( str_sample )
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
## All genes
lstr_write_genes = set( dict_gene_confirmed.keys() + dict_gene_missed.keys() )
## If key genes were given just look at key genes
if lstr_key_genes:
  lstr_write_genes = lstr_key_genes
## Counts of exome and RNA-Seq genes
li_counts = []
## Used to help make indices
i_cur_count_index = 0
## Indices for counts (bars)
li_counts_indices_exome = []
li_counts_indices_rna = []
## Count colors
lstr_plot_genes = []
li_exome_counts = []
li_rna_counts = []
## For each gene
for str_gene in lstr_write_genes:
  if dict_gene_confirmed.get( str_gene, 0 ) + dict_gene_missed.get( str_gene, 0 ):
    # In RNA-Seq
    li_rna_counts.append( dict_gene_confirmed.get( str_gene, 0 ) )
    # In exome
    li_exome_counts.append( dict_gene_confirmed.get( str_gene, 0 ) + dict_gene_missed.get( str_gene, 0 ) )
    # Update positioning
    i_cur_count_index = i_cur_count_index + 1
    li_counts_indices_rna.append( i_cur_count_index )
    i_cur_count_index = i_cur_count_index + 1
    li_counts_indices_exome.append( i_cur_count_index )
    i_cur_count_index = i_cur_count_index + 1
    # Add gene name
    lstr_plot_genes.append( str_gene )

# Plot counts for genes
## Positions for gene names
li_genes_indicies = [ ( i_index * 3 ) + 1.5 for i_index in range( len( lstr_plot_genes ) ) ]
fig, ax = plt.subplots()
ax.barh( li_counts_indices_rna, li_rna_counts, color=["cyan"]*len( li_counts_indices_rna), label="RNA-Seq" )
ax.barh( li_counts_indices_exome, li_exome_counts, color=["magenta"]*len( li_counts_indices_exome ), label ="Exome"  )
plt.yticks( li_genes_indicies, lstr_plot_genes )
ax.set_title( "SNPs in Exome Seq Recovered in RNA-Seq" )
ax.set_xlabel( "SNP (Count)" )
ax.set_ylabel( "Gene" )
for tick in ax.yaxis.get_major_ticks():
  tick.label.set_fontsize(6)
plt.legend( loc = 4 )
plt.savefig( args_call.str_output_file )
plt.close( fig )
## Write data to file
str_gene_text = os.path.splitext( args_call.str_output_file )[ 0 ] + ".txt"
with open( str_gene_text, "w" ) as hndl_gene:
  hndl_gene.write( "Gene\tExome\tRNA_Seq\n" )
  for i_gene_index in xrange( len( lstr_plot_genes ) ):
    hndl_gene.write( lstr_plot_genes[ i_gene_index ] + "\t" + str( li_exome_counts[ i_gene_index ] ) + "\t" + str( li_rna_counts[ i_gene_index ] ) + "\n" )

# Plot counts for genes
li_sample_counts_rna = []
li_sample_counts_exome = []
i_cur_count_index_sample = 0
li_count_indices_sample_rna = []
li_count_indices_sample_exome = []
lstr_plot_genes_sample = []
li_exome_counts_sample = []
li_rna_counts_sample = []
## For each gene
for str_gene in lstr_write_genes:
  if len( set( dict_gene_confirmed_sample.get( str_gene, [] ) + dict_gene_missed_sample.get( str_gene, [] ) ) ):
    # In RNA-Seq
    li_sample_counts_rna.append( len( set( dict_gene_confirmed_sample.get( str_gene, [] ) ) ) )
    li_rna_counts_sample.append( len( set( dict_gene_confirmed_sample.get( str_gene, [] ) ) ) )
    # In exome
    li_sample_counts_exome.append( len( set( dict_gene_confirmed_sample.get( str_gene, [] ) + dict_gene_missed_sample.get( str_gene, [] ) ) ) )
    li_exome_counts_sample.append( len( set( dict_gene_confirmed_sample.get( str_gene, [] ) + dict_gene_missed_sample.get( str_gene, [] ) ) ) )
    # Update positioning
    i_cur_count_index_sample = i_cur_count_index_sample + 1
    li_count_indices_sample_rna.append( i_cur_count_index_sample )
    i_cur_count_index_sample = i_cur_count_index_sample + 1
    li_count_indices_sample_exome.append( i_cur_count_index_sample )
    i_cur_count_index_sample = i_cur_count_index_sample + 1
    # Add gene name
    lstr_plot_genes_sample.append( str_gene )

# Plot counts for genes
## Positions for gene names
str_sample_plot = os.path.splitext( args_call.str_output_file )[ 0 ] + "_samples.pdf"
li_genes_indicies_samples = [ ( i_index * 3 ) + 1.5 for i_index in range( len( lstr_plot_genes_sample ) ) ]
fig, ax = plt.subplots()
ax.barh( li_count_indices_sample_rna, li_sample_counts_rna, color=["cyan"]*len(li_count_indices_sample_rna), label="RNA-Seq" )
ax.barh( li_count_indices_sample_exome, li_sample_counts_exome, color=["magenta"]*len(li_count_indices_sample_exome ), label="Exome" )
plt.yticks( li_genes_indicies_samples, lstr_plot_genes_sample )
ax.set_title( "SNPs in Exome Seq Recovered in RNA-Seq" )
ax.set_xlabel( "Sample (Count)" )
ax.set_ylabel( "Gene" )
for tick in ax.yaxis.get_major_ticks():
  tick.label.set_fontsize(6)
plt.legend( loc = 4 )
plt.savefig( str_sample_plot )
plt.close( fig )
## Write data to file
str_sample_text = os.path.splitext( str_sample_plot )[ 0 ] + ".txt"
with open( str_sample_text, "w" ) as hndl_sample:
  hndl_sample.write( "Gene\tExome_count\tRNA_Seq_count\n" )
  for i_sample_index in xrange( len( lstr_plot_genes_sample ) ):
    hndl_sample.write( lstr_plot_genes_sample[ i_sample_index ] + "\t" + str( li_exome_counts_sample[ i_sample_index ] ) + "\t" + str( li_rna_counts_sample[ i_sample_index ] ) + "\n" )
