#!/usr/bin/env python

# Import libraries
import argparse
import csv
import gzip
import os

# Constants
# VCF related
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"
I_INFO_INDEX = 7
CHR_INFO_DELIMITER = ";"
# INFO features
STR_CHASM_FDR = "CHASM_FDR"
STR_CHASM_PVALUE = "CHASM_PVALUE"
STR_VEST_FDR = "VEST_FDR"
STR_VEST_PVALUE = "VEST_PVALUE"
STR_COSMIC_ID = "COSMIC_ID"
STR_FATHMM = "FATHMM"
STR_CANCER = "CANCER"
# Thresholds
I_FDR = 0.3
I_PVALUE = 0.05

# Arguments
prog_desc = "".join(["Filters VCF based on mutation priority predictions."])
prsr_arguments = argparse.ArgumentParser(prog="filter_vcf_for_predictions.py",
                                         description=prog_desc,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
prsr_arguments.add_argument("str_input_file", help="Input vcf.gz file.")
prsr_arguments.add_argument("str_output_file", help="Output filtered vcf file.")
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []

def func_split_token(str_token):
    """
    Splits a VCF info token while guarenteeing 2 tokens
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

            # Do not keep any line that does not meet a criteria.
            f_keep = False

            # Keep comments
            if lstr_line[0][0] == CHR_COMMENT:
                lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
                continue

            # Parse INFO tokens
            dict_info_tokens = dict([func_split_token(str_token) for str_token
                                    in lstr_line[I_INFO_INDEX].split(CHR_INFO_DELIMITER)])

            # Keep everything that is a COSMIC ID at this point.
            # Otherwise require CRAVAT or VEST to have an annotation.
            if STR_CHASM_FDR in dict_info_tokens:
                if( float(dict_info_tokens.get(STR_CHASM_FDR, "2")) <= I_FDR):
                    f_keep = True
            elif(float(dict_info_tokens.get(STR_CHASM_PVALUE, "2")) <= I_PVALUE):
                f_keep = True
            if STR_VEST_FDR in dict_info_tokens:
                if( float(dict_info_tokens.get(STR_VEST_FDR, "2")) <= I_FDR):
                    f_keep = True
            elif(float(dict_info_tokens.get(STR_VEST_PVALUE, "2")) <= I_PVALUE):
                f_keep = True

            # Keep FATHMM = Cancer
            if( STR_COSMIC_ID in dict_info_tokens and
                dict_info_tokens.get(STR_FATHMM, "") == STR_CANCER ):
                f_keep = True

            # Store passing variant
            if f_keep:
                lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
                f_keep = False

        # Last write of buffer
        hndl_out.write("\n".join(lstr_vcf))

# Close input handle
hndl_vcf.close()
