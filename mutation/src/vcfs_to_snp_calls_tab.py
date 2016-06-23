#!/usr/bin/env python

# Constant
STR_POS_SEP = "--"
CHR_GENOTYPE_SEP = "/"

# VCF contants
I_CHR_INDEX = 0
I_POS_INDEX = 1
I_REF_INDEX = 3
I_ALT_INDEX = 4
I_FILTER_INDEX = 6
I_FORMAT_INDEX = 8
I_GENOTYPE_INDEX = 9
I_SINGLE_VCF_SAMPLE = 10
CHR_VCF_DELIM = "\t"
STR_PASS = "PASS"
STR_FORMAT_GENOTYPE = "GT"
CHR_ANNOT_DELIM = ":"
CHR_COMMENT = "#"
CHR_COMMENT_CHR = "#CHROM"
CHR_GENOTYPE_DELIM = "/"
CHR_GENOTYPE_DELIM_2 = "|"
CHR_GENOTYPE_MUT = ","
CHR_MONOMORPHIC_REFERENCE = "."

# Depth file constants
CHR_READ_DEPTH_DELIM = "\t"

# output file constants
CHR_OUTPUT_DELIM = "\t"
I_LOC_NOT_REF_VCF = 4
I_REF_NOT_REF_VCF = 5
I_GENOTYPE_NOT_REF_VCF = 6

# MAF constants
CHR_MAF_DELIM = "\t"
I_MAF_CHR_INDEX = 4
I_MAF_POS_INDEX = 5
I_MAF_REF_INDEX = 10
I_MAF_ALLELE_ONE = 11
I_MAF_ALLELE_TWO = 12
I_MAF_TUMOR_SAMPLE = 15
I_MAF_VARIANT_TYPE = 9
STR_MAF_SNP = "SNP"

import argparse
import csv
import glob
import gzip
import os

prsr_arguments = argparse.ArgumentParser( prog = "vcfs_to_snp_calls_tab", description = "Collapses vcfs to tabular files of snp calls", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-v", "--vcf", required = True, default = None, dest = "str_vcf_file", action = "store", help = "Input vcf file; maf or vcf are required." )
prsr_arguments.add_argument( "-r", "--vcf_reference", default = None, dest = "str_vcf_reference_file", action = "store", help = "Input reference vcf file; Any SNP in this file will be recorded even if it does not show up in the other VCF." )
prsr_arguments.add_argument( "-m", "--maf_reference", default = None, dest = "str_maf_file", action = "store", help = "Input reference maf file; reference maf or reference vcf are required." )
prsr_arguments.add_argument( "-t", "--tumor", default = None, dest = "str_maf_tumor_sample", action = "store", help = "The specific sample of interest to pull from the maf files, required when the reference file is a maf file.")
prsr_arguments.add_argument( "-d", "--count_reference", required = True, dest = "str_depth_reference_file", action = "store", help = "Input file that has per base read counts for the bam the reference vcf or maf was derived from.")
prsr_arguments.add_argument( "-c", "--count", required = True, dest = "str_vcf_depth_file", action = "store", help = "Input file that has per base read counts. For the bam used to drived the second vcf file which is being compared to the reference vcf / maf file.")
prsr_arguments.add_argument( "--no_prefilter", action="store_true", dest="f_no_prefilter_mode", help="When the VCF file has not been filtered (indicated by using this flag), this script will not require the SNPs to pass a filter")
prsr_arguments.add_argument( dest = "str_output_file", action = "store", help = "Output tab file of snp calls." )
args = prsr_arguments.parse_args()


def func_get_genotype( str_REF, str_ALT, str_genotype_info ):
  lstr_genotype = str_genotype_info.split( CHR_GENOTYPE_DELIM ) if CHR_GENOTYPE_DELIM in str_genotype_info else str_genotype_info.split( CHR_GENOTYPE_DELIM_2 )
  li_genotype = [ int( x ) for x in lstr_genotype ]
  if not len( li_genotype ) == 2:
    print "This script assumes diploid calls but received a genotype that is not diploid. entry = " + str_genotype_info
    exit( 11 )

  # Update
  # Sometimes the REF and ALT calls can be comma delimited
  lstr_REF = str_REF.split(",") if "," in str_REF else [ str_REF ]
  lstr_ALT = str_ALT.split(",") if "," in str_ALT else [ str_ALT ]
  if sum( [ 0 if len( str_ref_value ) == 1 else 1 for str_ref_value in lstr_REF ] ) > 0:
    return None
  lstr_genotype_values = lstr_REF + lstr_ALT

  # Make sure the reference position is 1 nucleotide only
  if not ( len( lstr_REF ) == 1 ):
    return None

  # Make sure the call indicates a change occured.
  if (( lstr_genotype_values[ li_genotype[0] ] == lstr_REF[ 0 ] ) and
     ( lstr_genotype_values[ li_genotype[1] ] == lstr_REF[ 0 ] )):
    return None

  # Make call
  str_genotype = CHR_GENOTYPE_SEP.join([ lstr_genotype_values[ li_genotype[0] ], lstr_genotype_values[ li_genotype[1] ]])

  # Check to make sure none of the calls are '.'
  if "." in str_genotype:
    return None

  # Make sure that the calls only have snps
  if not len( str_genotype ) == 3:
    return None
  return str_genotype


def func_update_with_depth( str_depth_file, dict_to_update, i_column_to_update ):
  print "Reading depth file: " + str_depth_file
  hndl_count = None
  if str_depth_file.split(".")[-1] == "gz":
    hndl_count = gzip.open( str_depth_file, "r" )
  else:
    hndl_count = open( str_depth_file, "r" )
  csvr_rd = csv.reader( hndl_count, delimiter = CHR_READ_DEPTH_DELIM )
  for lstr_count_line in csvr_rd:
    str_rd_loc = lstr_count_line[0]+"--"+lstr_count_line[1]
    if str_rd_loc in dict_output:
      dict_to_update[ str_rd_loc ][i_column_to_update] = lstr_count_line[2]
  hndl_count.close()
  return dict_to_update


def func_read_VCF_file( str_file_name, dict_reference = None, f_no_prefilter_mode = False ):

  dict_return = {}
  f_reference = not dict_reference == None
  if f_reference:
    dict_return = dict_reference

  print "Reading VCF file: " + str_file_name
  hndl_vcf = None
  if str_file_name.split( "." )[ -1 ] == ".gz":
    hndl_vcf = gzip.open( str_file_name, "r" )
  else:
    hndl_vcf = open( str_file_name, "r" )

  csvr_vcf = csv.reader( hndl_vcf, delimiter = CHR_VCF_DELIM )
  for lstr_line in csvr_vcf:
    # Skip comment
    if lstr_line[0][0] == CHR_COMMENT:
      if lstr_line[0] == CHR_COMMENT_CHR:
        if not len( lstr_line ) == I_SINGLE_VCF_SAMPLE:
          print "This VCF file does not adhere to the expected 1 sample format. Terminating."
      continue

    # Make sure the line passes
    if not f_no_prefilter_mode:
      if not lstr_line[I_FILTER_INDEX] in [ STR_PASS, "." ]:
        continue

    # Skip monomorphic sites
    if lstr_line[ I_ALT_INDEX ] == CHR_MONOMORPHIC_REFERENCE:
      continue

    # Get call information
    # Get the location of the call information from the format entry
    i_genotype_index = lstr_line[ I_FORMAT_INDEX ].split( CHR_ANNOT_DELIM ).index( STR_FORMAT_GENOTYPE )

    str_genotype = func_get_genotype( lstr_line[ I_REF_INDEX], lstr_line[ I_ALT_INDEX ], lstr_line[ I_GENOTYPE_INDEX ].split( CHR_ANNOT_DELIM )[ i_genotype_index ] )
    if not str_genotype:
      continue

    # Make genomic location
    str_loc = lstr_line[I_CHR_INDEX] + STR_POS_SEP + lstr_line[ I_POS_INDEX ]

    # Store in dict
    # f_reference indicates a reference dict has been given.
    # Both reference and not are recorded.
    if f_reference:
      if not str_loc in dict_return:
        dict_return[ str_loc ] = [ "NA" ] * 8
      dict_return[ str_loc ][ I_LOC_NOT_REF_VCF ] = str_loc
      dict_return[ str_loc ][ I_REF_NOT_REF_VCF ] = lstr_line[ I_REF_INDEX ]
      dict_return[ str_loc ][ I_GENOTYPE_NOT_REF_VCF ] = str_genotype
    else:
      dict_return[ str_loc ] = [str_loc, lstr_line[I_REF_INDEX], str_genotype, "NA", "NA", "NA", "NA", "NA" ]

  hndl_vcf.close()
  return dict_return

# Dict of output {chr--1232:line}
dict_output = None
# Read in the reference dict
if args.str_vcf_reference_file:

  # Read in the non-reference vcf file
  dict_output = func_read_VCF_file( str_file_name = args.str_vcf_reference_file, f_no_prefilter_mode=args.f_no_prefilter_mode )

elif args.str_maf_file:

  # The maf file is a consolidated file of results from many files. These files are using a potentially different naming convention.
  # This requires the name of the sample in the maf file of interest to be given for maf file parsing.
  if not args.str_maf_tumor_sample:
    print "Please supply the tumor sample of interest when parsing results from the maf file"
    exit( 12 )

  with open( args.str_maf_file, "r" ) as hndl_maf:

    print "Reading MAF file: " + args.str_maf_file
    print "Reporting on sample: " + args.str_maf_tumor_sample

    dict_output = {}

    csvr_maf = csv.reader( hndl_maf, delimiter = CHR_MAF_DELIM )
    for lstr_line_maf in csvr_maf:

      # Only look at targeted tumor sample
      if not lstr_line_maf[ I_MAF_TUMOR_SAMPLE ] == args.str_maf_tumor_sample:
        continue

      # Skip not snps
      if not lstr_line_maf[ I_MAF_VARIANT_TYPE ] == STR_MAF_SNP:
        continue

      # Get call
      # Make sure there is a difference in alleles and reference
      if (( lstr_line_maf[ I_MAF_REF_INDEX ] == lstr_line_maf[ I_MAF_ALLELE_ONE ]) and (lstr_line_maf[ I_MAF_ALLELE_ONE ] == lstr_line_maf[ I_MAF_ALLELE_TWO ])):
        continue
      str_call_maf = lstr_line_maf[ I_MAF_ALLELE_ONE ] + CHR_GENOTYPE_SEP + lstr_line_maf[ I_MAF_ALLELE_TWO ]

      # Make genomic location
      str_loc_maf = "" if lstr_line_maf[ I_MAF_CHR_INDEX ][0] == "c" else "chr"
      str_loc_maf = str_loc_maf + lstr_line_maf[ I_MAF_CHR_INDEX ] + STR_POS_SEP + lstr_line_maf[ I_MAF_POS_INDEX ]

      # Store in dict
      dict_output[ str_loc_maf ] = [ str_loc_maf, lstr_line_maf[ I_MAF_REF_INDEX ], str_call_maf, "NA", "NA", "NA", "NA", "NA" ]
else:
  print "Please provide either a reference VCF or maf file."
  exit( 13 )

# Read the not reference vcf
dict_output = func_read_VCF_file( str_file_name = args.str_vcf_file, dict_reference = dict_output, f_no_prefilter_mode=args.f_no_prefilter_mode )

# Dict to hold the read depth
dict_output = func_update_with_depth( str_depth_file = args.str_depth_reference_file, dict_to_update = dict_output, i_column_to_update = 3 )
dict_output = func_update_with_depth( str_depth_file = args.str_vcf_depth_file, dict_to_update = dict_output, i_column_to_update = 7 )

# Write to output file
with open( args.str_output_file, "w" ) as hndl_out:
  print "Writing output file: " + args.str_output_file
  # Write header
  hndl_out.write( "\n".join( [ CHR_READ_DEPTH_DELIM.join(lstr_values) for lstr_values in dict_output.values() ] ) )
