#!/usr/bin/env python

# Import libraries
import argparse
import csv
import gzip
import os

# Constants
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"
I_INFO_INDEX = 7
CHR_INFO_DELIMITER = ";"
STR_COSMIC_ID = "COSMIC_ID"
STR_COMMON_VARIANT = "COMMON"
STR_DEPTH = "DP"
STR_ORIGIN = "SAO"
STR_ORIGIN_GERMLINE = "1"
STR_FATHMM = "FATHMM"
STR_FATHMM_NEUTRAL = "PASSENGER/OTHER"

# Arguments
prog_desc = "".join(["Extracts common variants which do",
                     "not have COSMIC ids entries from a vcf file"])
prsr_arguments = argparse.ArgumentParser(prog="filter_vcf_for_cancer.py",
                                         description=prog_desc,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
prsr_arguments.add_argument("str_input_file", help="Input vcf.gz file.")
prsr_arguments.add_argument("str_output_file", help="Output filtered vcf file.")
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []
# Buffer size in lines
i_write_amount = 1000
# Count the filtering per filter
i_common = 0
i_depth = 0
i_origin = 0
i_fathmm = 0

def func_split_token(str_token):
  """
  Splits a VCF info token while guarenteeing 2 derived pieces.
  * str_token : String token to split at '='
              : String

  * return : List of 2 strings
  """
  if str_token:
    lstr_pieces = str_token.split("=")
    i_pieces = len(lstr_pieces)
    if i_pieces == 1:
      return(lstr_pieces + [""])
    if i_pieces == 2:
      return lstr_pieces
    elif i_pieces > 2:
      return [lstr_pieces[0], "=".join(lstr_pieces[1:])]
  return ["", ""]

# Get handle
str_input_ext = os.path.splitext(args.str_input_file)[1]
hndl_vcf = gzip.open(args.str_input_file, "rb") if str_input_ext == ".gz" else open(args.str_input_file, "rb")

# Read in vcf file
if args.str_input_file:
  with open(args.str_output_file, "w") as hndl_out:
    for lstr_line in csv.reader(hndl_vcf, delimiter = STR_VCF_DELIMITER):

      # Keep comments
      if lstr_line[0][0] == CHR_COMMENT:
        lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
        continue

      # Parse INFO tokens
      dict_info_tokens = dict([func_split_token(str_token) for str_token in lstr_line[I_INFO_INDEX].split(CHR_INFO_DELIMITER)])

      # Filtering
      ## Filter out common variant that do not have cosmic ids
      if STR_COMMON_VARIANT in dict_info_tokens:
        if((dict_info_tokens[STR_COMMON_VARIANT] == "1")
           and (STR_COSMIC_ID not in dict_info_tokens)):
          i_common = i_common + 1
          continue

      ## Filter DP < 10
      if STR_DEPTH in dict_info_tokens:
        if((int(dict_info_tokens[STR_DEPTH]) < 10)
           and (STR_COSMIC_ID not in dict_info_tokens)):
          i_depth = i_depth + 1
          continue

      ## Filter out SAO = Germline unless in COSMIC
      if STR_ORIGIN in dict_info_tokens:
        if((dict_info_tokens[STR_ORIGIN] == STR_ORIGIN_GERMLINE)
           and (STR_COSMIC_ID not in dict_info_tokens)):
          i_origin = i_origin + 1
          continue

      ## Filter out FATHMM = Neutral or passenger unless in COSMIC
      if STR_FATHMM in dict_info_tokens:
        if((dict_info_tokens[STR_FATHMM] == STR_FATHMM_NEUTRAL)
           and (STR_COSMIC_ID not in dict_info_tokens)):
          i_fathmm = i_fathmm + 1
          continue

      # Store passing variant
      lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))

      # If buffer is large enough, write to file
      if len(lstr_vcf) >= i_write_amount:
        str_write = "\n".join(lstr_vcf)
        if str_write:
          str_write = str_write + "\n"
        hndl_out.write(str_write)
        lstr_vcf = []

    # Last write of buffer
    hndl_out.write("\n".join(lstr_vcf))

# Close input handle
hndl_vcf.close()
