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

class MetagenomicsScript( PipelineRunner.PipelineRunner ):
    
    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "metagenomics.py "
        arg_raw.description = "Metagenomics - Centrifuge"
        arg_raw.add_argument("--threads", dest = "threads", default = "1", required = True ,help = "Launch NTHREADS parallel search threads - default: 1" )
        arg_raw.add_argument("--format",dest="format",choices=["fastq","fasta"],default="fastq" ,help="Choose format")
        arg_raw.add_argument("--num_primary_assign",dest="distinct_primary_assignments",default="1",help="It searches for at most <int> distinct, primary assignments for each read or pair.Default=5") 
        arg_raw.add_argument("--index", dest="centrifuge_index",help="The basename of the index for the reference genomes")
        arg_raw.add_argument("--read_type", dest="read_type", choices=["single","paired"], default="paired", help="Choose read type, skip if using Trinity assembles reads")
        arg_raw.add_argument("--right_fq", dest="right_fq" ,help="Right_fq (only when fastq format is used for read_type paired)")
        arg_raw.add_argument("--left_fq", dest="left_fq", help="Left_fq (only when fastq format is used for read_type paired)")
        arg_raw.add_argument("--unpaired_reads", dest="unpaired_reads",help="Comma-separated list of files containing unpaired reads to be aligned (for Trinity runs and single end reads)")
        return(arg_raw)

    def func_make_commands( self, args_parsed, cur_pipeline ):
        
        
        #For trinity assembled reads
        if args_parsed.format=="fasta":
            cur_pipeline.func_check_files_exist( [ args_parsed.unpaired_reads ] )
        #Direct route: paired or single end fastq 
        elif args_parsed.format=="fastq":
            if args_parsed.read_type=="single":
                cur_pipeline.func_check_files_exist( [ args_parsed.unpaired_reads ] )
            else:
                cur_pipeline.func_check_files_exist( [ args_parsed.right_fq,args_parsed.left_fq ] )
        else:
            cur_pipeline.logr_logger.error("Choose either --fasta or --fastq and provide appropriate inputs")

        # Make all output files 
        classification_results = os.path.join(args_parsed.str_out_dir,"classification.results.txt") 
        classification_report = os.path.join(args_parsed.str_out_dir,"classification.report.txt")
        report_kraken = os.path.join(args_parsed.str_out_dir,"kraken_style_report.txt") 
        
        #Trinity centrifuge command
        if args_parsed.format=="fasta":
            centrifige_cmd=["centrifuge",
                            "-p",args_parsed.threads,
                            "-k",args_parsed.distinct_primary_assignments,
                            "-f","-x",args_parsed.centrifuge_index,
                            "-U",args_parsed.unpaired_reads,
                            "-S",classification_results,
                            "--report-file",classification_report]

        #Direct fastq run with centrifuge
        elif args_parsed.format=="fastq":
            if args_parsed.read_type=="single":
                centrifige_cmd=["centrifuge",
                                "-p",args_parsed.threads,
                                "-k",args_parsed.distinct_primary_assignments,
                                "-q","-x",args_parsed.centrifuge_index,
                                "-U",args_parsed.unpaired_reads,
                                "-S",classification_results,
                                "--report-file",classification_report]
            else:
                centrifige_cmd=["centrifuge",              
                                 "-p",args_parsed.threads,
                                 "-q","-x",args_parsed.centrifuge_index,
                                 "-k",args_parsed.distinct_primary_assignments,
                                 "-1",args_parsed.left_fq,
                                 "-2",args_parsed.right_fq,
                                 "-S",classification_results,
                                 "--report-file",classification_report] 
        else:
            cur_pipeline.logr_logger.error("Choose either --fasta or --fastq and provide appropriate inputs")


        #Generate kraken style report
        kreport_cmd=["centrifuge-kreport",
                      "-x",args_parsed.centrifuge_index,
                      classification_results,">",report_kraken]
            
        
        lcmd_commands = []
        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join(centrifige_cmd),
                                                 lstr_cur_dependencies = [ args_parsed.centrifuge_index ], 
                                                 lstr_cur_products = [ classification_results,
                                                                       classification_report ] ) ] )

        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join(kreport_cmd),
                                                 lstr_cur_dependencies = [ args_parsed.centrifuge_index,classification_results],
                                                 lstr_cur_products = [ report_kraken ] ) ] )
        return(lcmd_commands)
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    MetagenomicsScript( ).func_run_pipeline( )
