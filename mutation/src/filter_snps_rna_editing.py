#!/usr/bin/env python

import argparse
import csv

# Constants
CHR_COMMENT="#"
STR_VCF_DELIMITER="\t"
I_CHR_INDEX=0
I_POS_INDEX=1
I_REF_INDEX=3
I_ALT_INDEX=4

STR_DARNED_DELIMITER="\t"
I_CHR_DARNED=0
I_POS_DARNED=1
I_STRAND_DARNED=2
I_INCHR_DARNED=3
I_INRNA_DARNED=4

STR_RADAR_DELIMITER="\t"
I_CHR_RADAR=0
I_POS_RADAR=1
I_STRAND_RADAR=3

str_prog_des = " ".join(["Filters SNP entries from a VCF file which are",
                         "in provided RNA Editing databases."])
fmtr_cur = argparse.ArgumentDefaultsHelpFormatter
prsr_arguments = argparse.ArgumentParser(prog="filter_snps_rna_editing.py",
                                         description=str_prog_des,
                                         formatter_class=fmtr_cur)
prsr_arguments.add_argument("--darned",
                            action="store",
                            dest="str_darned_db",
                            help=" ".join(["Darned database for RNA editing",
                                           "matched to the correct species",
                                           "reference."]))
prsr_arguments.add_argument("--radar",
                            action="store",
                            dest="str_radar_db",
                            help=" ".join(["Radar database for RNA editing",
                                           "matched to the correct species",
                                           "reference."]))
prsr_arguments.add_argument("str_input_file",
                            help="Input vcf file.")
prsr_arguments.add_argument("str_output_file",
                            help="Output SNP vcf file.")
args = prsr_arguments.parse_args()

# RNA editing
dict_darned = {}
dict_radar = {}
dict_seq_comp = {"A":"T","T":"A","G":"C","C":"G"}
dict_seq_interp = {"A":"A","C":"C","G":"G","T":"T","I":"G","U":"T"}

# Read in the RADAR data if given
if args.str_radar_db:
  with open(args.str_radar_db, "r") as hndl_radar:
    for lstr_line in csv.reader(hndl_radar, delimiter=STR_RADAR_DELIMITER):

      if not lstr_line:
        continue

      if lstr_line[I_STRAND_RADAR] == "+":
        str_rpos = lstr_line[I_CHR_RADAR].lower()+"-"+lstr_line[I_POS_RADAR]
        dict_radar[str_rpos] = ["A", dict_seq_interp["I"]]

      elif lstr_line[I_STRAND_RADAR] == "-":
        str_rneg = lstr_line[I_CHR_RADAR].lower()+"-"+lstr_line[I_POS_RADAR]
        dict_radar[str_rneg] = [dict_seq_comp["A"],
                                dict_seq_comp[dict_seq_interp["I"]]]

# Read in the DARNED data if given
if args.str_darned_db:
  with open(args.str_darned_db, "r") as hndl_darned:
    csvr_darned = csv.reader(hndl_darned, delimiter=STR_DARNED_DELIMITER)
    csvr_darned.next()
    for lstr_line in csvr_darned:

      if not lstr_line:
        continue

      inchr = lstr_line[I_INCHR_DARNED]
      inrna = lstr_line[I_INRNA_DARNED]

      if lstr_line[I_STRAND_DARNED] == "+":
        str_dpos = lstr_line[I_CHR_DARNED].lower()+"-"+lstr_line[I_POS_DARNED]
        dict_darned[str_dpos] = [dict_seq_interp[inchr],
                                 dict_seq_interp[inrna]]

      elif lstr_line[I_STRAND_DARNED] == "-":
        str_dmin = lstr_line[I_CHR_DARNED].lower()+"-"+lstr_line[I_POS_DARNED]
        dict_darned[str_dmin] = [dict_seq_comp[dict_seq_interp[inchr]],
                                 dict_seq_comp[dict_seq_interp[inrna]]]

# Stores the vcf info
lstr_vcf = []

# Read in vcf file
if args.str_input_file:
  with open(args.str_output_file, "w") as hndl_out:
    with open(args.str_input_file, "r") as hndl_vcf:
      for lstr_line in csv.reader(hndl_vcf, delimiter=STR_VCF_DELIMITER):

        # Keep comments
        if lstr_line[0][0] == CHR_COMMENT:
          lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
          continue

        # Get ALT / REF
        str_alt = lstr_line[I_ALT_INDEX]
        str_ref = lstr_line[I_REF_INDEX]

        # Create ID (chr-pos)
        str_id = lstr_line[I_CHR_INDEX].lower()+"-"+lstr_line[I_POS_INDEX]

        # Filter Darned
        if dict_darned:
          if str_id in dict_darned:
            lstr_rna_edit_entry = dict_darned[str_id]
            if(lstr_rna_edit_entry[0] == str_ref
               and lstr_rna_edit_entry[1] == str_alt):
              continue

        # Filter Radar
        if dict_radar:
          if str_id in dict_radar:
            lstr_rna_edit_entry = dict_radar[str_id]
            if(lstr_rna_edit_entry[0] == str_ref
               and lstr_rna_edit_entry[1] == str_alt):
              continue

        # Store SNP
        lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))

      for str_out_line in lstr_vcf:
        hndl_out.write(str_out_line + "\n")
