#!/usr/bin/env python


__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2014"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"


import argparse
import rnaseq_mutation_pipeline as rnaseq


# Parse arguments
prsr_arguments = argparse.ArgumentParser( prog = "rnaseq_mutation_indexer.py", description = "Creates initial star aligner index for the rnaseq mutation pipeline", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
prsr_arguments.add_argument( "-b", "--bsub_queue", metavar = "BSUB_queue", dest = "str_bsub_queue", default = None, help = "If given, each command will sequentially be ran on this queue with bsub." )
prsr_arguments.add_argument( "-c", "--clean", dest = "f_clean", default = False, action="store_true", help = "Turns on (true) or off (false) cleaning of intermediary product files." ) 
prsr_arguments.add_argument( "-d", "--alignment_mode", metavar = "Alignment_mode", dest = "str_alignment_mode", default = rnaseq.STR_ALIGN_STAR, choices = rnaseq.LSTR_ALIGN_CHOICES, help = "Specifies the alignment and indexing algorithm to use." )
prsr_arguments.add_argument( "-f", "--reference", metavar = "Reference genome", dest = "str_genome_fa", required = True, help = "Path to the reference genome to use in the analysis pipeline." )
prsr_arguments.add_argument( "-g", "--log", metavar = "Optional_logging_file", dest = "str_log_file", default = None, help = "Optional log file, if not given logging will be to the standard out." )
prsr_arguments.add_argument( "-k", "--gtf", metavar = "Reference GTF", dest = "str_gtf_file_path", default = None, help = "GTF file for reference genome.")
prsr_arguments.add_argument( "-m", "--max_bsub_memory", metavar = "Max_BSUB_memory", dest = "str_max_memory", default = "8", help = "The max amount of memory requested in GB when running bsub commands." )
prsr_arguments.add_argument( "-n", "--threads", metavar = "Process_threads", dest = "i_number_threads", type = int, default = 1, help = "The number of threads to use for multi-threaded steps." )
prsr_arguments.add_argument( "-o", "--out_dir", metavar = "Output_directory", dest = "str_file_base", default = None, help = "The output directory where results will be placed. If not given a directory will be created from sample names and placed with the samples." )
prsr_arguments.add_argument( "-t", "--test", dest = "f_Test", default = False, action = "store_true", help = "Will check the environment and display commands line but not run.")
args = prsr_arguments.parse_args()

# Needed for compatibility with pipeline
# If the output file is not given it is based off
# Normally the name of the left file but in this case
# the name of the reference genome.
args.str_sample_file_left_fq = args.str_genome_fa

rnaseq.run( args, True )
