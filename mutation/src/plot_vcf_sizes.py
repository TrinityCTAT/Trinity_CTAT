#!/usr/bin/env python

# Import libraries
import argparse
import csv
import os
import barChart
import quickPlot as qp

# Constants
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"

# Arguments
prsr_arguments = argparse.ArgumentParser( prog = "plot_vcf_sizes.py", description = "Plot the number of features per vcf file", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-v","--vcf", dest="lstr_input_vcfs", action="append", help = "Input vcf files to measure ( can be a mixture of vcf )." )
prsr_arguments.add_argument( "-o", "--output", dest="str_pdf", action="store", help = "Output plot pdf." )
args = prsr_arguments.parse_args()

# Hold the number of of noncommment features per file
lstr_labels = []
li_counts = []

# Get handle
for str_vcf in args.lstr_input_vcfs:

  i_feature_count = 0
  with open( str_vcf, "rb" ) as hndl_vcf:
    for lstr_line in csv.reader( hndl_vcf, delimiter = STR_VCF_DELIMITER ):

      # Keep comments
      if lstr_line[0][0] == CHR_COMMENT:
        continue

      # Count line
      i_feature_count = i_feature_count + 1

  # Store feature info
  lstr_labels.append( os.path.basename( str_vcf ) )
  li_counts.append( i_feature_count )

# Plot feature count of each file.
barChart.BarChart().func_plot( { qp.c_STR_TITLE : "Impact of Filtering per Step",
                     qp.c_STR_X_AXIS : "Filter Step",
                     qp.c_STR_Y_AXIS : "Number of Features Left in File",
                     qp.c_STR_DATA : [ { qp.c_STR_DATA : li_counts,
                                         qp.c_STR_X_TICK_LABEL : lstr_labels } ] }, args.str_pdf )
