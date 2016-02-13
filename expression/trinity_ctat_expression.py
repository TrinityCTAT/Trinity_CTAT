#!/usr/bin/env python


__author__ = "Asma Bankapur"
__copyright__ = "Copyright 2014"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "bankapur@broadinstitute.org"
__status__ = "Development"

import os
import sciedpiper.Command as Command
import sciedpiper.ParentScript as ParentScript

KALLISTO = 'kallisto'
kallisto_script = 'kallisto_script'

class ExpressionScript( ParentScript.ParentScript ):
    
    
    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.description = "Expression from RNASeq Data for Trinity CTAT"
        arg_raw.add_argument( "--left_fq", required=True, help = "Left read fastq" )        
        arg_raw.add_argument( "--right_fq", default="", help ="Right read fastq" )
        arg_raw.add_argument( "--annot_config", required=True ,help ="Annotation Config file" )
        arg_raw.add_argument( "--bias", action="store_true", help= "Perform sequence based bias correction")
        arg_raw.add_argument( "--bootstrap_samples", help ="Number of bootstrap samples" )
        arg_raw.add_argument( "--threads", help ="Number of threads to use for bootstraping" )
        arg_raw.add_argument( "--seed", help ="Seed for the bootstrap sampling" )
        

    def func_make_commands( self, args_parsed, cur_pipeline ):
        
        # Make directories and check files that need to exist before beginning
 
        kallisto_dir = os.path.join( args_parsed.str_file_base, KALLISTO )
        cur_pipeline.func_mkdirs( [ kallisto_dir  ] )

        ### Get kallisto script command
        kallisto_cmd = [ kallisto_script + ' quant' ] 
        
        if args_parsed.annot_config is not None:
           with open( args_parsed.annot_config, 'r' ) as annot:
                for line in annot:
                    annot_list = line.split("=") 
                    if annot_list[ 0 ] == 'KALLISTO':
                       kallisto_index = annot_list[ 1 ] 
           option_cmd1 = [ " --index " + kallisto_index.strip("\n") ]
           kallisto_cmd.extend( option_cmd1 )
           
        if args_parsed.bias:
           option_cmd2 = [ " --bias " ]
           kallisto_cmd.extend( option_cmd2 )

        if not args_parsed.bootstrap_samples:
           option_cmd3 = [ " --bootstrap-samples " + args_parsed.bootstrap_samples ]
           kallisto_cmd.extend( option_cmd3 )

        if not args_parsed.threads:
           option_cmd4 = [ " --threads " + args_parsed.threads ]
           kallisto_cmd.extend( option_cmd4 )

        if not args_parsed.seed:
           option_cmd5 = [ " --seed  " + args_parsed.seed ]
           kallisto_cmd.extend( option_cmd5 )
        
        option_cmd6 = [ " -o " + kallisto_dir ]
        kallisto_cmd.extend( option_cmd6 )

        if not args_parsed.right_fq:
           option_cmd7 = [ " --single " + args_parsed.left_fq ] 
        else:
           option_cmd7 = [ args_parsed.left_fq + " " + args_parsed.right_fq ]
      
        kallisto_cmd.extend( option_cmd7 )
        kallisto_cmd_str = " ".join( kallisto_cmd )        


        #print kallisto_cmd_str        
        ###Final list of commands to run 
        lcmd_commands = [ ]
        lcmd_commands.extend( [ Command.Command( str_cur_command = kallisto_cmd_str,
                                                 lstr_cur_dependencies = [ args_parsed.left_fq ], 
                                                 lstr_cur_products = [ kallisto_dir ] ) ] )
        return lcmd_commands
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    ExpressionScript().func_run_pipeline()
