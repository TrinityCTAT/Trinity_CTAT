#!/usr/bin/env python

# Constants
CHR_COMMENT = "#"
CHR_COMMENT_CHR = "#CHROM"
STR_VCF_DELIMITER = "\t"
STR_COMMON = "COMMON"
I_CHR_INDEX = 0
I_SINGLE_VCF_SAMPLE = 10
I_POS_INDEX = 1
I_REF_INDEX = 3
I_ALT_INDEX = 4
I_FORMAT_INDEX = 7
I_GENOTYPE_INDEX = 9
I_FILTER_INDEX = 6
STR_PASS = "PASS"
CHR_MONOMORPHIC_REFERENCE = "."
STR_FORMAT_GENOTYPE = "GT"
STR_SOMATIC = "SAO"
CHR_ANNOT_DELIM = ";"
STR_POS_SEP = "--"

import argparse
import csv
import os

prsr_arguments = argparse.ArgumentParser( prog = "annotate_vcf.py", description = "Annotates a vcf file given different sources", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-d","--dbsnp", dest = "str_dbsnp", help = "DBSNP vcf file of annotation." )
prsr_arguments.add_argument( "-i", "--input", dest = "lstr_input_files", action = "append", help = "Input vcf file." )
prsr_arguments.add_argument( "str_output_file", help = "Output annotated vcf file." )
args = prsr_arguments.parse_args()

# Annotate with dbsnp
lstr_comments = []
# All the mutations in the reference vcf file
# ( assumed to be prefiltered )
# { chr--pos: {REF: {ALT : [ "ANNOTATION,ANNOTATION,ANNOTATION" ]} } }
dict_mutations = {}
if args.str_dbsnp:
  with open( args.str_input_file, "r" ) as hndl_dbsnp_vcf:
    for lstr_line in csv.reader( hndl_dbsnp_vcf, delimiter = STR_VCF_DELIMITER ):
      
      # Store comment in comment list
      if lstr_line[0][0] == CHR_COMMENT:
        lstr_comments.append( lstr_line )
        continue

      # Store genomic feature
      str_pos = lstr_line[ I_CHR_INDEX ] + STR_POS_SEP + lstr_line[ I_POS_INDEX ]

      # Annotation info
      str_feature_annotation = lstr_line[ I_FORMAT_INDEX ]
      dict_feature_annotation_tokens = dict( [ str_annotation_token.split("=") for str_annotation_token in str_feature_annotation.split(";") ] )

      # Add COMMON annotation
      lstr_bases = lstr_line[ I_REF ].split(",") + lstr_line[ I_ALT ].split(",")
      # If common
      if( "COMMON=1" in str_feature_annotation ):
        # If known percent
        if( "CAF=[" in str_feature_annotation ):
          lstr_caf = dict_feature_annotation_tokens[ "CAF" ].strip("[]").split(",")
          for i_index in range( 1 : len( lstr_bases )):
            dict_mutations[ str_pos ] = { lstr_bases[ 0 ]: { lstr_bases[ i_index ]:[ "COMMON(=" + lstr_caf[ i_index ] + ")" ] } }
        # If unknown percent
        else
          for str_alt in lstr_bases[1:]:
            dict_mutations[ str_pos ] = { lstr_bases[ 0 ]: { str_alt: [ "COMMON(>=0.01)" ] } }
        # If not common
        elif( "COMMON=0" in str_feature_annotation ):
          for str_alt in lstr_bases[1:]:
            dict_mutations[ str_pos ] = { lstr_bases[ 0 ]: { str_alt: [ "NOT_COMMON(<0.01)" ] } }
        else:
          for str_alt in lstr_bases[1:]:
            dict_mutations[ str_pos ] = { lstr_bases[ 0 ]: { str_alt: [ "COMMON_UNKNOWN" ] } }
      
       # SAO = 0
       if( "SAO" in str_feature_annotation ):
         str_somatic_key = "SOMATIC_STATUS:"
         # SAO = 0 Unspecified
         if( "SAO=0" in str_feature_annotation ):
           str_somatic_key = str_somatic_key + "UNSPECIFIED"
         # SAO = 1
         elif( "SAO=1" in str_feature_annotation ):
           str_somatic_key = str_somatic_key + ""
         # SAO = 2
         elif( "SAO=2" in str_feature_annotation ):
           str_somatic_key = str+somatic_key + ""
         # SAO = 3 BOTH
         elif( "SAO=3" in str_feature_annotation ):
           str_somatic_key = str_somatic_key + "BOTH"
         else:
           str_somatic_key = str_somatic_key + "UNSPECIFIED"

        # Maybe I will just keep all the DBSNP info
        # Indicate if there is a 1000 Genome submission for it KGPROD
        # Indicate if it is validated in 1000 Genome KGVALIDATED
        # Cited ? MUT
        # PMC ? Pubmed
        # PM precious
        # DBSNP id ( rs number )
        # WTD?
        # VLD? for mult alleles

        # Indicate somatic mutations
        # Seems that is
        ##INFO=<ID=SAO,Number=1,Type=Integer,Description=""Variant Allele Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both"">
        if lstr_token[ 0 ] == STR_SOMATIC:
          str_variant_type = lstr_token[ 1 ]

        # Store Allelic frequency from 1000Genomes CAF

        # Store Variant Disease name CLNDBN

        # CLNSIG variant clinical significance

        # Gene symbols GENEINFO

      # Make genomic location
      str_cur_loc = lstr_line[I_CHR_INDEX] + STR_POS_SEP + lstr_line[ I_POS_INDEX ]

      if str_prev_loc == str_cur_loc:
        print "INFO::" + str_cur_loc + " has more than one entry."
      str_prev_loc = str_cur_loc
      llstr_list.append( str_cur_loc + ".." + str_ref +  "->" + str_alt + "\t" + str_variant_type + "\n" )

# Read in the VCF file
# and write vcf file
#with open( str_input_file, "r" ) as hndl_in:
#  # Write vcf file
#  with open( str_output_file, "w" ) as hndl_out:
#    for str_line in hndl_in:
#      # Write header
#      hndl_out.write( str_line )
