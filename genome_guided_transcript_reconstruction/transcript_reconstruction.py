#!/usr/bin/env python


__author__ = "Asma Bankapur, Timothy Tickle"
__copyright__ = "Copyright 2014"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

#import inspect
import os, sys
sys.path.append(os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "SciEDPipeR"]))
sys.path.append(os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "SciEDPipeR", "sciedpiper"]))
import sciedpiper.Command as Command
import sciedpiper.PipelineRunner as PipelineRunner

OUTPUT_GTF_MINUS_COMMENTS = "transcripts_reconstructed_minus_comments.gtf"

class TransReconstruction( PipelineRunner.PipelineRunner ):
    
    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "transcript_reconstruction.py"
        arg_raw.description = "Assembling Transcripts - STRINGTIE"
        arg_raw.add_argument( "--bam_file", required=True ,help ="Aligned Bam file" )
        arg_raw.add_argument( "--ref_annot", required=True , help="Reference annotation, GTF" )
        return(arg_raw) 

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Allows:
        - the creation of commands in the child object.
        - the creation of directories.
        - checking that files exist.
        
        To know the variables available from command line look in the PipelineRunner in func_create_arguments.
        """
        
        cur_pipeline.func_check_files_exist( [ args_parsed.bam_file,args_parsed.ref_annot ] )
       
        output_gtf=os.path.join(args_parsed.str_out_dir, "transcripts.gtf")

        output_bed=os.path.join(args_parsed.str_out_dir, "transcripts.bed")
 
        output_intermediate=os.path.join(args_parsed.str_out_dir, OUTPUT_GTF_MINUS_COMMENTS)
        # Make strigtie command:
        stringtie_cmd_list = [ 'stringtie', 
                               args_parsed.bam_file,
                               '-G', args_parsed.ref_annot,
                               '-o', output_gtf ]
        
       
        stringtie_cmd = " ".join( stringtie_cmd_list )

        # Strip comments from stringtie's gtf

        strip_cmd_list = [ 'tail -n +3',
                           output_gtf, '>',
                           output_intermediate ]

        strip_cmd = " ".join( strip_cmd_list ) 


        # Make gtf2bed command:
        gtf_to_bed_cmd_list  = [ 'gtf2bed.py', 
                                 output_intermediate,
                                 '>', output_bed ]

        gtf_to_bed_cmd = " ".join( gtf_to_bed_cmd_list ) 

        lcmd_commands = []
        lcmd_commands.extend( [ Command.Command( str_cur_command = stringtie_cmd,
                                                 lstr_cur_dependencies = [ args_parsed.bam_file, args_parsed.ref_annot ], 
                                                 lstr_cur_products = [ output_gtf ] ),
                                
                                Command.Command( str_cur_command = strip_cmd,
                                                 lstr_cur_dependencies = [ output_gtf ],
                                                 lstr_cur_products = [ output_intermediate ] ),

                                Command.Command( str_cur_command = gtf_to_bed_cmd,
                                                 lstr_cur_dependencies = [ output_intermediate ] ,
                                                 lstr_cur_products = [ output_bed ] )
                              ] ) 
                                                                     
        return(lcmd_commands)
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    TransReconstruction( ).func_run_pipeline( )
