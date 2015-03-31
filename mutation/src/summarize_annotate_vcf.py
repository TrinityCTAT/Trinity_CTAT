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

class SummarizeAnnotateVCF( ParentScript.ParentScript ):

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
    args_raw.add_argument( "--output_file", metavar="output_file", dest="str_output_file", required=True, action="store", help="Final output vcf file (please include vcf as the extension." )
    args_raw.add_argument( metavar = "sample_vcfs", dest = "lstr_sample_files", default = None, action="store", nargs="+",  help = "Sample files to combine vcfs." )
    print( vars(args_raw.parse_args()) )

  def func_make_commands( self, args_parsed, cur_pipeline ):
    """
    Make a list of command objects.

    * return : List of commands to be ran
             : List of command objects
    """

    # File of combined sample vcf files
    str_combined_vcf = os.path.splitext( args_parsed.str_output_file )[ 0 ] + "_combined.vcf"

    # Manage Files
    # Make sure the output file is named gz, if not, add it
    if not os.path.splitext( args_parsed.str_output_file )[ 1 ] == "gz":
      args_parsed.str_output_file = args_parsed.str_output_file + ".gz"

    # List of commands
    lcmd_commands = []

    # Combine sample vcf files if more than one sample vcf file is given
    # bcftools merge --merge all --output_type z --output str_output_file.vcf.gz lstr_vcf_files
    if len( args_parsed.lstr_sample_files ) > 2:
        str_merge_command = " ".join( [ "bcftools", "merge", "--merge", "all", "--output-type", "z", "--output" + str_combined_vcf ] + args_parsed.lstr_sample_files )
        lcmd_commands.append( Command.Command( str_cur_command = str_merge_command,
                                               lstr_cur_dependencies = args_parsed.lstr_sample_files,
                                               lstr_cur_products = [ str_combined_vcf ] ) )
    else:
        str_combined_vcf = args_parsed.lstr_sample_files[ 0 ]
        # If using an input file compress and index
        if not os.path.splitext( str_combined_vcf )[ 1 ] == ".gz":
            lcmd_commands.append( Command.Command( str_cur_command = "bgzip -c "+str_combined_vcf + " > " + str_combined_vcf + ".gz",
                                                   lstr_cur_dependencies = [ str_combined_vcf ],
                                                   lstr_cur_products = [ str_combined_vcf+".gz" ] ) )
            str_combined_vcf = str_combined_vcf + ".gz"
            lcmd_commands.append( Command.Command( str_cur_command = "tabix -p vcf " + str_combined_vcf,
                                                   lstr_cur_dependencies = [ str_combined_vcf ],
                                                   lstr_cur_products = [ str_combined_vcf+".tbi" ] ) )

    if not os.path.splitext( args_parsed.str_dbsnp_vcf )[ 1 ] == ".gz":
        lcmd_commands.append( Command.Command( str_cur_command = "bgzip -c " + args_parsed.str_dbsnp_vcf + " > " + args_parsed.str_dbsnp_vcf + ".gz",
                                                   lstr_cur_dependencies = [ args_parsed.str_dbsnp_vcf ],
                                                   lstr_cur_products = [ args_parsed.str_dbsnp_vcf + ".gz" ] ) )
        args_parsed.str_dbsnp_vcf = args_parsed.str_dbsnp_vcf + ".gz"
        lcmd_commands.append( Command.Command( str_cur_command = "tabix -p vcf " + args_parsed.str_dbsnp_vcf,
                                                   lstr_cur_dependencies = [ args_parsed.str_dbsnp_vcf ],
                                                   lstr_cur_products = [ args_parsed.str_dbsnp_vcf + ".tbi" ] ) )
    # Annotate combined sample vcf files
    # bcftools annotate --annotations str_dbsnp_vcf -c
    # PM variant is clinicall precious (clinical and pubmed cited)
    # NSF, NSM, NSN, COMMON, SAO, KGPROD, KGVALIDATED, MUT, WTD, VLD, RS, PMC
    str_annotate_command = " ".join( [ "bcftools", "annotate", "--output-type", "z", "--annotations", args_parsed.str_dbsnp_vcf, "--columns", "INFO/COMMON,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SAO,INFO/KGPROD,INFO/KGValidated,INFO/MUT,INFO/WTD,INFO/VLD,INFO/RS,INFO/PMC", "--header-lines", args_parsed.str_add_headers_file, "--output", args_parsed.str_output_file, str_combined_vcf ] )
    lcmd_commands.append( Command.Command( str_cur_command = str_annotate_command,
                                           lstr_cur_dependencies = [ args_parsed.str_dbsnp_vcf, str_combined_vcf ],
                                           lstr_cur_products = [ args_parsed.str_output_file ] ) )

    # Run commands
    return( lcmd_commands )

if __name__ == "__main__":

    # Needed to run, calls the script
    SummarizeAnnotateVCF().func_run_pipeline()
