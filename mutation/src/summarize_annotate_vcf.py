#!/usr/bin/env python

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2015"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

import argparse
import os
import sciedpiper.Command as Command
import sciedpiper.Pipeline as Pipeline
import sciedpiper.ParentScript as ParentScript

class SummarizeAnnotateVCF( ParentScript ):

  def func_update_arguments( self, args_raw ): 
    """
    Add in custom arguments fo rthe mutation validation pipeline.

    * arg_raw : Arguments ( not yet parsed )
              : Arguments
    * return  : Updated Arguments
              : Arguments
    """

    # Parse command line arguments
    args_raw.prog = "Summarize_annotate_vcf.py"
    args_raw.add_argument( "--dbsnp", metavar = "dbsnp_reference_vcf", dest = "str_dbsnp_vcf", default = None, required = True, action="store", help = "Reference dbsnp file to use for annotation." )
    args_raw.add_argument( "--reference_headers", metavar = "add_headers_file", dest = "str_add_headers_file", required = True, action="store", help = "A test file of headers to add (info being added in annotations needs a header." )
    args_raw.add_argument( metavar = "sample_vcfs", dest = "lstr_sample_files", default = None, required = True, action="store", nargs="+",  help = "Sample files to combine vcfs." )
    return( args_raw )


  def func_make_commands( self, args_parsed, cur_pipeline ):
    """
    Make a list of command objects.

    * return : List of commands to be ran
             : List of command objects
    """

    # File of combined sample vcf files
    str_combined_vcf = os.path.splitext( args_call.str_output_file )[0] + "_combined.vcf"

    # Manage Files
    # Make sure the output file is named gz, if not, add it
    if not os.path.splitext( args_call.str_output_file )[ 1 ] == "gz":
      args_call.str_output_file = args_call.str_output_file + ".gz"

    # List of commands
    lcmd_commands = []

    # Combine sample vcf files
    # bcftools merge --merge all --output_type z --output str_output_file.vcf.gz lstr_vcf_files
    str_merge_command = " ".join( [ "bcftools", "merge", "--merge", "all", "--output_type", "z", "--output" + str_output_file.vcf.gz ] + lstr_vcf_files )
    lcmd_commands.append( Command.Command( str_cur_command = str_merge_command,
                     lstr_cur_dependencies = lstr_vcf_files,
                     lstr_cur_products = [ str_combined_vcf ] ) )

    # Annotate combined sample vcf files
    # bcftools annotate --annotations str_dbsnp_vcf -c
    # PM variant is clinicall precious (clinical and pubmed cited)
    # NSF, NSM, NSN, COMMON, SAO, KGPROD, KGVALIDATED, MUT, WTD, VLD, RS, PMC
    str_annotate_command = " ".join( [ "bcftools", "annotate", "--annotations", str_dbsnp_vcf, "--columns", "INFO/COMMON,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SAO,INFO/KGPROD,INFO/KGVALIDATED,INFO/MUT,INFO/WTD.INFO/VLD,INFO/RS,INFO/PMC", "--header_lines", args_call.str_add_headers_file, args_call.str_output_file ] )
    lcmd_commands.append( Command.Command( str_cur_command = str_annotate_command,
                     lstr_cur_dependencies = [ str_dbsnp_vcf, str_combined_vcf ],
                     lstr_cur_products = [ args_call.str_output_file ] ) )

    # Run commands
    return( lcmd_commands )
