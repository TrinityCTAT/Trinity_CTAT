#!/usr/bin/env python

# Constants
CHR_COMMENT = "#"
C_STR_ANN_HEADER = "##INFO=<ID=ANN"
C_STR_ANN_HEADER_PREFIX = "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: '"
STR_VCF_DELIMITER = "\t"
C_GENE_NAME_SNPEFF = "Gene_Name"
C_GENE_NAME = "GENE"
C_STR_CHROM = "#CHROM"
C_STR_GENE_NAME_HEADER = "##INFO=<ID="+C_GENE_NAME+",Number=.,Type=String,Description=\"The name of the gene/s in the genomic region of the SNP as annotated by SNPeff\">"
C_STR_ANN = "ANN"
C_I_INFO_INDEX = 7

import argparse
import csv

prsr_arguments = argparse.ArgumentParser( prog = "update_snpeff_annotations.py", description = "Updated SNPEff annotations.", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "str_input_file", help = "Input vcf file." )
prsr_arguments.add_argument( "str_output_file", help = "Output vcf file with updated annotations." )
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []

# Index fo the Gene_Name in the SNPEff annotation
i_gene_name = -1

# Read in vcf file
if args.str_input_file:
  with open( args.str_input_file, "r" ) as hndl_vcf:
    for lstr_line in csv.reader( hndl_vcf, delimiter = STR_VCF_DELIMITER ):

      # Store comments
      if lstr_line[0][0] == CHR_COMMENT:
        str_current_comment_line = STR_VCF_DELIMITER.join( lstr_line )

        # If the comment is the SNPEff ANN header element, find the location of the Gene_Name
        if str_current_comment_line.split(",")[0] == C_STR_ANN_HEADER:
          lstr_ann_tokens = [ str_token.strip() for str_token in str_current_comment_line[ len( C_STR_ANN_HEADER_PREFIX ) : ].split("|") ]
          i_gene_name = lstr_ann_tokens.index( C_GENE_NAME_SNPEFF )

        # Add the new Gene_name info line
        if str_current_comment_line.split(STR_VCF_DELIMITER)[0] == C_STR_CHROM:
         lstr_vcf.append( C_STR_GENE_NAME_HEADER  )

        # Store comment
        lstr_vcf.append( str_current_comment_line )
        continue

      # Parse for the potentially multiple ANN entries and gene names contained within
      lstr_ann = [ str_token for str_token in lstr_line[ C_I_INFO_INDEX ].split(";") if str_token.split("=")[0] == C_STR_ANN ]
      str_ann_gene_names = "NA"
      if len( lstr_ann ) == 1:
        str_ann = lstr_ann[0]
        lstr_ann_gene_names = [ str_token.split("|")[ i_gene_name ] for str_token in str_ann.split(",") ]
        str_ann_gene_names = ",".join( set(lstr_ann_gene_names))

      # Add gene names to info block
      lstr_line[ C_I_INFO_INDEX ] = C_GENE_NAME+"="+str_ann_gene_names+";"+lstr_line[ C_I_INFO_INDEX ]

      # Store body
      lstr_vcf.append( STR_VCF_DELIMITER.join( lstr_line ) )
      continue

with open( args.str_output_file, "w" ) as hndl_out:
  for str_out_line in lstr_vcf:
      hndl_out.write( str_out_line + "\n" )
