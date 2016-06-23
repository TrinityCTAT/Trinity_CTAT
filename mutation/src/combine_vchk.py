#!/usr/bin/env python

import argparse
import barChart as bc
import boxPlot as bx
import glob
import json
import matplotlib.pyplot as plt
import os
import quickPlot as qp


def func_write_json( dict_json, str_file_name ):
    """ Write a dict representing json to a file """

    with open( str_file_name, "w" ) as hndl_out:
        hndl_out.write( json.dumps( dict_json, sort_keys=True, indent=2 ) )


prsr_arguments = argparse.ArgumentParser( prog = "combine_vchk.py", description = "Combines vchk files from samtools and makes some charts on combined statistics (Currently just the transitions data)", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "--input_dir", required = True, dest = "str_input_dir", action = "store", help = "Directory that contains vchk files (recursive glob will grab any file with .vchk in root or child dirs." )
prsr_arguments.add_argument( "--output_dir", required = True, dest = "str_output_dir", action = "store", help = "Output folder for text and plot files." )
args_call = prsr_arguments.parse_args()

# Check for the input directory
if not os.path.exists( args_call.str_input_dir ):
    print( "Error. The input dir does not exist. " + args_call.str_input_dir )
    exit( 1 )

# Get the .vchk files
lstr_vchk_files = []
for str_root, lstr_dir, lstr_files in os.walk( args_call.str_input_dir ):
    lstr_vchk_files.extend( [ os.path.join( str_root, str_file ) for str_file in lstr_files if os.path.splitext( str_file )[ 1 ] == ".vchk" ] )

# Stop if no files found
if len( lstr_vchk_files ) < 1:
    print( "No .vchk files found in the following directory: " + args_call.str_input_dir )
    exit( 2 )

# Make the output dir
if not os.path.exists( args_call.str_output_dir ):
    os.mkdir( args_call.str_output_dir )

# Parse the vchk files
dict_absolute = {}
dict_percent = {}
for str_vchk_file in lstr_vchk_files:
    with open( str_vchk_file ) as hndl_vchk:
        print( str_vchk_file )
        f_read = False
        dict_file_values = {}
        i_file_total = 0.0
        for str_line in hndl_vchk:
            if not f_read:
                if str_line == "# ST, Substitution types:\n":
                    f_read = True
            else:
                if not str_line[ 0 ] == "#":
                    lstr_current_tokens = [ str_token for str_token in str_line.split("\t") if str_token ]
                    if( not len( lstr_current_tokens ) == 4 ):
                        print( "Error parsing this line in this file" )
                        print( str_line )
                        print( str_vchk_file )
                        hndl_vchk.close()
                        exit( 3 )
                    else:
                        dict_file_values[ lstr_current_tokens[ 2 ] ] = int( lstr_current_tokens[ 3 ] )
                        i_file_total = i_file_total + int( lstr_current_tokens[ 3 ] )
        # Update global dict totals
        for str_key, i_value in dict_file_values.items():
            dict_percent.setdefault( str_key, [] ).append( i_value / i_file_total )
            dict_absolute[ str_key ] = dict_absolute.setdefault( str_key, 0 ) + i_value

# Set up the dicts for plotting
lstr_labels = []
lli_values = []
for str_label, livalues in dict_percent.items():
    lstr_labels.append( str_label )
    lli_values.append( livalues )

# Set up dict for absolute
lstr_labels_abs = dict_absolute.keys()
lli_values_abs = [ dict_absolute[ str_key ] for str_key in lstr_labels_abs ]

# Plot
# Should move this to quickplots
# Boxplots
dict_rel_json = { qp.c_STR_TITLE: "Mutations by Substitution",
                  qp.c_STR_X_AXIS: "Base Substitution",
                  qp.c_STR_Y_AXIS: "Percent",
                  qp.c_STR_DATA: lli_values,
                  qp.c_STR_DATA_LABEL: lstr_labels }
bx.BoxPlot().func_plot( dict_rel_json, os.path.join( args_call.str_output_dir, "Distributions_substitutions.pdf" ))
func_write_json( dict_json=dict_rel_json, str_file_name= os.path.join( args_call.str_output_dir, "Distributions_substitutions.json" ) )

# Barchart
dict_absolute_json = {  qp.c_STR_TITLE: "Mutations by Substition",
                        qp.c_STR_X_AXIS: "Base Substitutions",
                        qp.c_STR_Y_AXIS: "Total Count",
                        qp.c_STR_DATA: [ { qp.c_STR_DATA : lli_values_abs,
                                           qp.c_C_PLOT_COLOR : qp.c_C_PLOT_COLOR_DEFAULT,
                                           qp.c_STR_X_TICK_LABEL : lstr_labels_abs,  } ] }
bc.BarChart().func_plot( dict_absolute_json, os.path.join( args_call.str_output_dir, "Total_substitutions.pdf" ))
func_write_json( dict_json=dict_absolute_json, str_file_name= os.path.join( args_call.str_output_dir, "Total_substitutions.json" ) )
