#!/usr/bin/env python


__author__ = "Asma Bankapur, Timothy Tickle"
__copyright__ = "Copyright 2014"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Asma Bankapur"
__email__ = "bankapur@broadinstitute.org"
__status__ = "Development"

#import inspect
import os,sys
sys.path.append(os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "SciEDPipeR"]))
sys.path.append(os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "SciEDPipeR", "sciedpiper"]))
import sciedpiper.Command as Command
import sciedpiper.PipelineRunner as PipelineRunner

class LncrnaScript( PipelineRunner.PipelineRunner ):
    
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
        arg_raw.add_argument("--config", type=str, help="path to assembly.config file. default uses config file in same directory as slncky")
        arg_raw.add_argument("--no_orth_search",action="store_true", help="flag if you only want to filter lncs but don\'t want to search for orthologs")
        arg_raw.add_argument('--no_filter', action='store_true', help='flag if you don\'t want lncs to be filtered before searching for ortholog')
        arg_raw.add_argument('--overwrite', action='store_true', help='forces overwrite of out_prefix.bed')
        arg_raw.add_argument('--threads', type=int, help='number of threads. default = 5', default=5)
        arg_raw.add_argument('--min_overlap', type=float, help='remove any transcript that overlap annotated coding gene > min_overlap%%. default = 0%%', default=0)
        arg_raw.add_argument('--min_cluster', type=int, help='min size of duplication clusters to remove. default=2', default=2)
        arg_raw.add_argument('--min_coding', type=float, help='min exonic identity to filter out transcript that aligns to orthologous coding gene. default is set by learning coding alignment distribution from data', default=0.1)
        arg_raw.add_argument('--no_overlap', action='store_true', help='flag if you don\'t want to overlap with coding')
        arg_raw.add_argument('--no_collapse', action='store_true', help='flag if you don\'t want to collapse isoforms')
        arg_raw.add_argument('--no_dup', action='store_true', help='flag if don\'t want to align to duplicates')
        arg_raw.add_argument('--no_self', action='store_true', help='flag if you don\'t want to self-align for duplicates')
        arg_raw.add_argument('--no_coding', action='store_true', help='flag if you don\'t want to align to orthologous coding')
        arg_raw.add_argument('--no_bg', action='store_true', help='flag if you don\'t want to compare lnc-to-ortholog alignments to a background. This flag may be useful if you want to do a \'quick-and-dirty\' run of the ortholog search.')
        arg_raw.add_argument('--no_orf', action='store_true', help='flag if you don\'t want to search for orfs')
        arg_raw.add_argument('--minMatch', type=float, help='minMatch parameter for liftover. default=0.1', default=0.1)
        arg_raw.add_argument('--pad', type=int, help='# of basepairs to search up- and down-stream when lifting over lnc to ortholog', default=50000)
        arg_raw.add_argument('--gap_open', type=str, default='200')
        arg_raw.add_argument('--gap_extend', type=str, default='40')
        arg_raw.add_argument('--web', action='store_true', help='flag if you want website written visualizing transcripts that were filtered out')        
        arg_raw.add_argument('--lastz', help='lastz software')
        arg_raw.add_argument('--bedtools', help='bedtools software')
        arg_raw.add_argument('--liftover',  help='liftover software')
        return(arg_raw)

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Allows:
        - the creation of commands in the child object.
        - the creation of directories.
        - checking that files exist.
        
        To know the variables available from command line look in the PipelineRunner in func_create_arguments.
        """

        
        # Bed file and config file must be present
        cur_pipeline.func_check_files_exist( [ args_parsed.bedfile,args_parsed.config ] )
        

        # Make web dir if --web included
        if args_parsed.web:
           web_dir = ( args_parsed.out_prefix + '.EvolutionBrowser' ) 
           cur_pipeline.func_mkdirs( [ web_dir ] )
   
        # Make all other files 
        # TODO: include conditionals based on usage
        canonical_to_lncs = args_parsed.out_prefix + ".canonical_to_lncs.txt" 
        cluster_info = args_parsed.out_prefix + ".cluster_info.txt"
        filtered_info = args_parsed.out_prefix + ".filtered_info.txt" 
        lncs_bed = args_parsed.out_prefix + ".lncs.bed" 
        lncs_info = args_parsed.out_prefix + ".lncs.info.txt" 
        orfs = args_parsed.out_prefix + ".orfs.txt" 
        orthologs_top = args_parsed.out_prefix + ".orthologs.top.txt" 
        orthologs = args_parsed.out_prefix + ".orthologs.txt" 
        
        # Make slncky command:
        slncky_cmd_list = [ 'slncky.v1.0', 
                            '--config', args_parsed.config,
                            '--threads', str( args_parsed.threads ),
                            '--min_overlap', str( args_parsed.min_overlap ),
                            '--min_cluster', str( args_parsed.min_cluster ),
                            '--min_coding', str( args_parsed.min_coding ),
                            '--bedtools',args_parsed.bedtools,
                            '--liftover',args_parsed.liftover,
                            '--minMatch',str( args_parsed.minMatch ),
                            '--pad',str( args_parsed.pad ),
                            '--lastz',args_parsed.lastz,
                            '--gap_open',str( args_parsed.gap_open ),
                            '--gap_extend',str( args_parsed.gap_extend ) ]
        
        boolean_args_list_all = [ '--no_orth_search','--no_filter',
                                  '--overwrite','--no_overlap',
                                  '--no_collapse','--no_dup',
                                  '--no_self','--no_coding',
                                  '--no_bg','--no_orf','--web' ]

        args_list_all =         [ args_parsed.no_orth_search,args_parsed.no_filter,
                                  args_parsed.overwrite,args_parsed.no_overlap,
                                  args_parsed.no_collapse,args_parsed.no_dup,
                                  args_parsed.no_self,args_parsed.no_coding,
                                  args_parsed.no_bg,args_parsed.no_orf,args_parsed.web ]

       
        for i,argument in enumerate( boolean_args_list_all ):
            if args_list_all[ i ] :
               slncky_cmd_list.append( argument )
        
        slncky_cmd_list.append( args_parsed.bedfile )
        slncky_cmd_list.append( args_parsed.assembly )
        slncky_cmd_list.append( args_parsed.out_prefix )        
       
        if slncky_cmd_list:
           slncky_cmd = " ".join( slncky_cmd_list )

        lcmd_commands = []
        lcmd_commands.extend( [ Command.Command( str_cur_command = slncky_cmd,
                                                 lstr_cur_dependencies = [ args_parsed.config ], 
                                                 lstr_cur_products = [ canonical_to_lncs ] ) ] )
        
        return(lcmd_commands)
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    LncrnaScript( ).func_run_pipeline( )
