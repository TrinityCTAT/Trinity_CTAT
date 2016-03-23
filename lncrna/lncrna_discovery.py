#!/usr/bin/env python


__author__ = "Asma Bankapur, Timothy Tickle"
__copyright__ = "Copyright 2014"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

#import inspect
import os
import sciedpiper.Command as Command
import sciedpiper.ParentScript as ParentScript

class LncrnaScript( ParentScript.ParentScript ):
    
    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "lncrna_discovery.py"
        arg_raw.description = "Lncrna and ortholog discovery"
        arg_raw.add_argument("--bedfile", dest = "bedfile", required = True ,help = "bed12 file of transcripts" )
        arg_raw.add_argument("--assembly" ,required = True ,help="assembly")
        arg_raw.add_argument("--out_prefix", default="slncky",help="out_prefix")
        arg_raw.add_argument("--config", default=CONFIGSTR ,type=str, help="path to assembly.config file. default uses config file in same directory as slncky")
        arg_raw.add_argument("--no_orth_search",action="store_true", help="flag if you only want to filter lncs but don\'t want to search for orthologs")
        arg_raw.add_argument('--no_filter', action='store_true', help='flag if you don\'t want lncs to be filtered before searching for ortholog')
        arg_raw.add_argument('--overwrite', action='store_true', help='forces overwrite of out_prefix.bed')
        arg_raw.add_argument('--threads', type=int, help='number of threads. default = 5', default=5)
        arg_raw.add_argument('--min_overlap', type=float, help='remove any transcript that overlap annotated coding gene > min_overlap%%. default = 0%%', default=0)
        arg_raw.add_argument('--min_cluster', type=int, help='min size of duplication clusters to remove. default=2', default=2)
        arg_raw.add_argument('--min_coding', type=float, help='min exonic identity to filter out transcript that aligns to orthologous coding gene. default is set by learning coding alignment distribution from data', default=None)
        arg_raw.add_argument('--no_overlap', action='store_true', help='flag if you don\'t want to overlap with coding')
        arg_raw.add_argument('--no_collapse', action='store_true', help='flag if you don\'t want to collapse isoforms')
        arg_raw.add_argument('--no_dup', action='store_true', help='flag if don\'t want to align to duplicates')
        arg_raw.add_argument('--no_self', action='store_true', help='flag if you don\'t want to self-align for duplicates')
        arg_raw.add_argument('--no_coding', action='store_true', help='flag if you don\'t want to align to orthologous coding')
        arg_raw.add_argument('--min_noncoding', type=float, help='min exonic identity to filter out transcript that aligns to orthologous noncoding gene. default=0', default=0.0)
        arg_raw.add_argument('--no_bg', action='store_true', help='flag if you don\'t want to compare lnc-to-ortholog alignments to a background. This flag may be useful if you want to do a \'quick-and-dirty\' run of the ortholog search.')
        arg_raw.add_argument('--no_orf', action='store_true', help='flag if you don\'t want to search for orfs')
        arg_raw.add_argument('--minMatch', type=float, help='minMatch parameter for liftover. default=0.1', default=0.1)
        arg_raw.add_argument('--pad', type=int, help='# of basepairs to search up- and down-stream when lifting over lnc to ortholog', default=0)
        arg_raw.add_argument('--gap_open', type=str, default='200')
        arg_raw.add_argument('--gap_extend', type=str, default='40')
        arg_raw.add_argument('--web', action='store_true', help='flag if you want website written visualizing transcripts that were filtered out')        


    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Allows:
        - the creation of commands in the child object.
        - the creation of directories.
        - checking that files exist.
        
        To know the variables available from command line look in the ParentScript in func_create_arguments.
        """

        
        # Make output directory
        cur_pipeline.func_mkdirs( [ args_parsed.str_file_base ] )
        # Bed file and config file must be present
        cur_pipeline.func_check_files_exist( [ args_parsed.bedfile,args_parsed.configstr ] )
        

        # Make web dir if --web included
        if args_parsed.web:
           web_dir = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + '.EvolutionBrowser' ) )
           cur_pipeline.func_mkdirs( [ web_dir ] )
           os.path.join( args_parsed.str_file_base, web_dir, "browse.html" )
   
        # Make all other files 
        # TODO: include conditionals based on usage
        canonical_to_lncs = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".canonical_to_lncs.txt" ) )
        cluster_info = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".cluster_info.txt" ) )
        filtered_info = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".filtered_info.txt" ) )
        lncs_bed = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".lncs.bed" ) )
        lncs_info = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".lncs.info.txt" ) )
        orfs = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".orfs.txt" ) )
        orthologs_top = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".orthologs.top.txt" ) )
        orthologs = os.path.join( args_parsed.str_file_base, ( args_parsed.out_prefix + ".orthologs.txt" ) )
        
        # Make slncky command:
        slncky_cmd_list = [ 'slncky.v1.0', 
                            '--config', args_parsed.config,
                            '--threads', args_parsed.threads,
                            '--min_overlap', args_parsed.min_overlap,
                            '--min_cluster', args_parsed.min_cluster,
                            '--min_coding', args_parsed.min_noncoding,
                            '--min_noncoding', args_parsed.min_noncoding,
                            '--bedtools',args_parsed.bedtools,
                            '--liftover',args_parsed.liftover,
                            '--minMatch',args_parsed.minMatch,
                            '--pad',args_parsed.pad,
                            '--lastz',args_parsed.lastz,
                            '--gap_open',args_parsed.gap_open,
                            '--gap_extend',args_parsed.gap_extend ]
        
        boolean_args_list_all = [ 'no_orth_search','no_filter',
                                  'overwrite','no_overlap',
                                  'no_collapse','no_dup',
                                  'no_self','no_coding',
                                  'no_bg','no_orf','web' ]

        args_dict = vars( args_parsed ) 
       
        for argument in boolean_args_list_all:
            if argument in args_dict:
                 slncky_cmd_list.append( "--" + argument )
        
        slncky_cmd_list.append( args_parsed.bedfile )
        slncky_cmd_list.append( args_parsed.assembly )
        slncky_cmd_list.append( args_parsed.out_prefix )        
       
        slncky_cmd = " ".join( slncky_cmd_list )

        lcmd_commands = []
        lcmd_commands.extend( [ Command.Command( str_cur_command = slncky_cmd,
                                                 lstr_cur_dependencies = [ args_parsed.config ], 
                                                 lstr_cur_products = [ canonical_to_lncs, cluster_info,
                                                                       filtered_info, lncs_bed,
                                                                       lncs_info, orfs,
                                                                       orthologs_top, orthologs ] ) ] )
        
        return lcmd_commands
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    LncrnaScript( ).func_run_pipeline( )
