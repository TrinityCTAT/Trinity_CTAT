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
    args_raw.add_argument( "--output_file", metavar="output_file", dest="str_output_file", required=True, action="store", help="Final output vcf file (please include vcf as the extension." )
    args_raw.add_argument( "--darned", metavar = "Darned_data", dest = "str_darned", default = None, action="store", help = "The DARNED database preped for annotation to be used in annotating RNA editing." )
    args_raw.add_argument( "--radar", metavar = "Radar_data", dest = "str_radar", default = None, action="store", help = "The RADAR database preped for annotation to beused in annotating RNA editing." )
    args_raw.add_argument( metavar = "sample_vcfs", dest = "lstr_sample_files", default = None, action="store", nargs="+",  help = "Sample files to combine vcfs." )

  def func_make_commands( self, args_parsed, cur_pipeline ):
    """
    Make a list of command objects.

    * return : List of commands to be ran
             : List of command objects
    """

    # File of combined sample vcf files
    str_combined_vcf = os.path.splitext( args_parsed.str_output_file )[ 0 ] + "_combined.vcf.gz"
    str_current_vcf = str_combined_vcf
    str_dbsnp_annotated_vcf = os.path.splitext( args_parsed.str_output_file )[ 0 ] + "_combined_dbsnp.vcf.gz"

    # List of commands
    lcmd_commands = []

    # List of files (compressed vcf) to work with
    lstr_compressed = []

    # If not compressed this will be updated
    str_annotation_file = args_parsed.str_dbsnp_vcf

    # Project directory
    str_project_dir = os.path.dirname( args_parsed.str_output_file )

    # Convert to compressed form if not
    for str_input_vcf in args_parsed.lstr_sample_files:
        if not os.path.splitext( str_input_vcf )[ 1 ] == ".gz":
            str_compressed_file = os.path.join( str_project_dir, os.path.basename( str_input_vcf ) + ".gz" )
            lcmd_commands.append( Command.Command( str_cur_command = "bgzip -c "+str_input_vcf + " > " + str_compressed_file,
                                               lstr_cur_dependencies = [ str_input_vcf ],
                                               lstr_cur_products = [ str_compressed_file ] ) )
            lcmd_commands.append( Command.Command( str_cur_command = "tabix -p vcf " + str_compressed_file,
                                               lstr_cur_dependencies = [ str_compressed_file ],
                                               lstr_cur_products = [ str_compressed_file + ".tbi" ] ) )
            lstr_compressed.append( str_compressed_file )
        else:
            lstr_compressed.append( str_input_vcf )

    # Convert dbsnp reference to compressed if needed
    if not os.path.splitext( str_annotation_file )[ 1 ] == ".gz":
        str_annotation_file = os.path.join( str_project_dir, os.path.basename( args_parsed.str_dbsnp_vcf ) + ".gz" )
        lcmd_commands.append( Command.Command( str_cur_command = "bgzip -c " + args_parsed.str_dbsnp_vcf + " > " + str_annotation_file,
                                                   lstr_cur_dependencies = [ args_parsed.str_dbsnp_vcf ],
                                                   lstr_cur_products = [ str_annotation_file ] ) )
        lcmd_commands.append( Command.Command( str_cur_command = "tabix -p vcf " + str_annotation_file,
                                                   lstr_cur_dependencies = [ str_annotation_file ],
                                                   lstr_cur_products = [ str_annotation_file + ".tbi" ] ) )

    # Plot vcf file
    for str_compressed in lstr_compressed:
        str_vchk_stats = os.path.splitext( str_compressed )[0] + ".vchk"
        str_plot_location = os.path.join( os.path.dirname( str_vchk_stats ), os.path.basename( str_vchk_stats ) + "_plot" )
        str_vchk_stats_command = " ".join( [ "bcftools", "stats", str_compressed, ">", str_vchk_stats ] )
        str_vchk_plot_command = " ".join( [ "plot-vcfstats", str_vchk_stats, "-p", str_plot_location + os.path.sep ] )
        lcmd_commands.append( Command.Command( str_cur_command = str_vchk_stats_command,
                                               lstr_cur_dependencies = [ str_compressed ],
                                               lstr_cur_products = [ str_vchk_stats ] ) )
        lcmd_commands.append( Command.Command( str_cur_command = str_vchk_plot_command,
                                               lstr_cur_dependencies = [ str_vchk_stats ],
                                               lstr_cur_products = [ str_plot_location ] ) )

    # Combine sample vcf files if more than one sample vcf file is given
    # bcftools merge --merge all --output_type z --output str_output_file.vcf.gz lstr_vcf_files
    if len( lstr_compressed ) > 1:
        str_merge_command = " ".join( [ "bcftools", "merge", "--merge", "all", "--output-type", "z", "--output", str_combined_vcf ] + lstr_compressed )
        lcmd_commands.append( Command.Command( str_cur_command = str_merge_command,
                                               lstr_cur_dependencies = lstr_compressed,
                                               lstr_cur_products = [ str_combined_vcf ] ) )
        lcmd_commands.append( Command.Command( str_cur_command = "tabix -p vcf " + str_combined_vcf,
                                                   lstr_cur_dependencies = [ str_combined_vcf ],
                                                   lstr_cur_products = [ str_combined_vcf + ".tbi" ] ) )
    else:
        str_combined_vcf = lstr_compressed[ 0 ]

    # DBSNP annotation
    # Annotate combined sample vcf files
    # bcftools annotate --annotations str_dbsnp_vcf -c
    # PM variant is clinicall precious (clinical and pubmed cited)
    # NSF, NSM, NSN, COMMON, SAO, KGPROD, KGVALIDATED, MUT, WTD, VLD, RS, PMC
    str_annotate_command = " ".join( [ "bcftools", "annotate", "--output-type", "v", "--annotations", str_annotation_file, "--columns", "INFO/COMMON,INFO/PM,INFO/NSF,INFO/NSM,INFO/NSN,INFO/SAO,INFO/KGPROD,INFO/KGValidated,INFO/MUT,INFO/WTD,INFO/VLD,INFO/RS,INFO/PMC", "--output", str_dbsnp_annotated_vcf, str_combined_vcf ] )
    lcmd_commands.append( Command.Command( str_cur_command = str_annotate_command,
                                           lstr_cur_dependencies = [ str_annotation_file, str_combined_vcf ],
                                           lstr_cur_products = [ str_dbsnp_annotated_vcf ] ) )
    str_current_vcf = str_dbsnp_annotated_vcf

#    # DARNED annotation
#    if args_parsed.str_darned:
#        str_darned_annotated_vcf = os.path.splitext( str_current_vcf )[0] + "_darned.vcf.gz" if args_parsed.str_radar else args_parsed.str_output_file
#        str_annotate_darned = " ".join( [ "bcftools","annotate", "--output-type", "v", "--columns", "CHR,POS,-,DARNED,-,-,-,-,-,-,-,-", "--annotations", args_parsed.str_darned, "--output", str_darned_annotated_vcf, str_current_vcf ] )
#        cmd_darned = Command.Command( str_cur_command = str_annotate_darned,
#                         lstr_cur_dependencies = [ str_current_vcf ],
#                         lstr_cur_products = [ str_darned_annotated_vcf ] )
#        lcmd_commands.append( cmd_darned )
#        str_current_vcf = str_darned_annotated_vcf
#
#    # RADAR annotation
#    if args_parsed.str_radar:
#        str_annotate_radar = " ".join( [ "bcftools", "annotate", "--output-type", "v", "--columns", "CHR,POS,RADAR,-,-,-,-,-,-,-,-", "--annotations", args_parsed.str_radar, "--output", args_parsed.str_output_file, str_current_vcf ] )
#        cmd_annotate_radar = Command.Command( str_cur_command = str_annotate_radar,
#                         lstr_cur_dependencies = [ str_current_vcf ],
#                         lstr_cur_products = [ args_parsed.str_output_file ] ) 
#        str_current_vcf = args_parsed.str_output_file

    # Run commands
    return( lcmd_commands )

if __name__ == "__main__":

    # Needed to run, calls the script
    SummarizeAnnotateVCF().func_run_pipeline()
