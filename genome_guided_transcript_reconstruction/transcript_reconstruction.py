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

OUTPUT_GTF_MINUS_COMMENTS = "transcripts_reconstructed_minus_comments.gtf"

class TransReconstruction( ParentScript.ParentScript ):
    
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
        arg_raw.add_argument( "--ref_annot", required=True , help="Reference annotation" )
        arg_raw.add_argument( "--output_gtf", required=True , help="Output GTF" )        
        arg_raw.add_argument( "--output_bed", required=True , help="Output GTF" )

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Allows:
        - the creation of commands in the child object.
        - the creation of directories.
        - checking that files exist.
        
        To know the variables available from command line look in the ParentScript in func_create_arguments.
        """
        
        cur_pipeline.func_check_files_exist( [ args_parsed.bam_file,args_parsed.ref_annot ] )
        
        # Make strigtie command:
        stringtie_cmd_list = [ 'stringtie', 
                               args_parsed.bam_file,
                               '-G', args_parsed.ref_annot,
                               '-o', args_parsed.output_gtf ]
        
       
        stringtie_cmd = " ".join( stringtie_cmd_list )

        # Strip comments from stringtie's gtf

        strip_cmd_list = [ 'tail -n +3',
                           args_parsed.output_gtf, '>',
                           OUTPUT_GTF_MINUS_COMMENTS ]

        strip_cmd = " ".join( strip_cmd_list ) 


        # Make gtf2bed command:
        gtf_to_bed_cmd_list  = [ 'gtf2bed.py', 
                                 OUTPUT_GTF_MINUS_COMMENTS,
                                 '>', args_parsed.output_bed ]

        gtf_to_bed_cmd = " ".join( gtf_to_bed_cmd_list ) 

        lcmd_commands = []
        lcmd_commands.extend( [ Command.Command( str_cur_command = stringtie_cmd,
                                                 lstr_cur_dependencies = [ args_parsed.bam_file, args_parsed.ref_annot ], 
                                                 lstr_cur_products = [ args_parsed.output_gtf ] ),
                                
                                Command.Command( str_cur_command = strip_cmd,
                                                 lstr_cur_dependencies = [ args_parsed.output_gtf ],
                                                 lstr_cur_products = [ OUTPUT_GTF_MINUS_COMMENTS ] ),

                                Command.Command( str_cur_command = gtf_to_bed_cmd,
                                                 lstr_cur_dependencies = [ OUTPUT_GTF_MINUS_COMMENTS ] ,
                                                 lstr_cur_products = [ args_parsed.output_bed ] )
                              ] ) 
                                                                     
        return lcmd_commands
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    TransReconstruction( ).func_run_pipeline( )
