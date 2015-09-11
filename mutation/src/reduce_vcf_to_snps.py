#!/usr/bin/env python

# Constants
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"
I_SINGLE_VCF_SAMPLE = 10
I_REF_INDEX = 3
I_ALT_INDEX = 4
I_FILTER_INDEX = 6
I_INFO_INDEX = 7
STR_PASS = "pass"
STR_VC_SNP = "snv"
CHR_MONOMORPHIC_REFERENCE = "."
CHR_INFO_DELIMITER = ";"

import argparse
import csv

prsr_arguments = argparse.ArgumentParser( prog = "reduce_vcf_to_snps.py", description = "Extracts only snp entries from a vcf file", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--reference", action = "store_true", dest = "f_reference_mode", help = "Reference vcf file mode. Does not check to see if the vcf feature passes." )
prsr_arguments.add_argument( "--no_prefilter", action="store_true", dest="f_no_prefilter_mode", help="When the VCF file has not been filtered (indicated by using this flag), this script will not require the SNPs to pass a filter")
prsr_arguments.add_argument( "str_input_file", help = "Input vcf file." )
prsr_arguments.add_argument( "str_output_file", help = "Output SNP vcf file." )
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []

i_write_amount = 1000

# Read in vcf file
if args.str_input_file:
  with open( args.str_output_file, "w" ) as hndl_out:
    with open( args.str_input_file, "r" ) as hndl_vcf:
      for lstr_line in csv.reader( hndl_vcf, delimiter = STR_VCF_DELIMITER ):

        # Keep comments
        if lstr_line[0][0] == CHR_COMMENT:
          lstr_vcf.append( STR_VCF_DELIMITER.join( lstr_line ) )
          continue
        # Look for reference vcf indication of mutation type
        if args.f_reference_mode:
          lstr_info = lstr_line[ I_INFO_INDEX ].split( CHR_INFO_DELIMITER )
          f_found_mutation = False
          for str_info in lstr_info:
            lstr_info_token = str_info.split("=")
            if lstr_info_token[ 0 ].lower() == "vc" and lstr_info_token[ 1 ].lower() == STR_VC_SNP:
              f_found_mutation = True
              continue
          if f_found_mutation:
            lstr_vcf.append( STR_VCF_DELIMITER.join( lstr_line ) )
            #print "Found SNP info tag, stored."
            continue

        # If is not a reference file, there should be a call for passing or not.
        # Make sure the variant passes.
        elif not args.f_no_prefilter_mode:
          if not lstr_line[ I_FILTER_INDEX ].lower() == STR_PASS.lower():
            #print "Did not pass filter. "+lstr_line[ I_FILTER_INDEX ].lower()+"."
            continue

        # Get ALT / REF
        str_alt = lstr_line[ I_ALT_INDEX ]
        str_ref = lstr_line[ I_REF_INDEX ]

        # Skip monomorphic sites
        if str_alt == CHR_MONOMORPHIC_REFERENCE or str_ref == CHR_MONOMORPHIC_REFERENCE:
          #print "Skip monomorph"
          continue

        # Skip if not SNPs
        # Want to keep anything that could be a SNP so
        # if the reference and the alt both have 1 base entries they should be kept, even if there are none SNP entries
        # Not removing the none snp information currently because the DBSNP info would have to be updated to reflect any indexing
        # referring to the REF or ALT
        if "," in str_alt:
          if not min( [ len( str_alt_token ) for str_alt_token in str_alt.split( "," ) ] ) == 1:
            #print "Skip not snp ALT with comma"
            continue
        elif len( str_alt ) > 1:
            #print "Skip not snp ALT"
            continue
        if "," in str_ref:
          if not min( [ len( str_ref_token ) for str_ref_token in str_ref.split( "," ) ] ) == 1:
            #print "Skip not snp REF with comma"
            continue
        elif len( str_ref ) > 1:
          #print "Skip not snp REF"
          continue

        # Store SNP
        lstr_vcf.append( STR_VCF_DELIMITER.join( lstr_line ) )

        # If buffer is large enough, write to file
        if len( lstr_vcf ) >= i_write_amount:
          for str_out_line in lstr_vcf:
            hndl_out.write( str_out_line + "\n" )
          lstr_vcf = []

      for str_out_line in lstr_vcf:
        hndl_out.write( str_out_line + "\n" )
