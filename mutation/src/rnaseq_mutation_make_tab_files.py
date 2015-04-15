#!/usr/bin/env python


__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2015"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

import csv
import os
import sciedpiper.Command as Command
import sciedpiper.Pipeline as Pipeline
import sciedpiper.ParentScript as ParentScript

class RNASEQ_mutation_make_tabs( ParentScript.ParentScript ):

    def __init__( self ):

        # File structure
        self.str_depth_dir = "depth"
        self.str_log_dir = "log"
        self.str_src_dir = "src"
        self.str_tab_dir = "tab"
        self.str_filtered_vcf = "vcf_filtered"

    def func_update_arguments( self, arg_raw ):
        """
        Add in custom arguments for the mutation validation pipeline.
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "RNASEQ_mutation_validation"

        # File that pairs RNA and DNA seq files
        # DNA BAM
        arg_raw.add_argument( "--dna_bam", dest = "str_dna_bam", required=True, help = "Input DNA BAM file." )
        # DNA VCF
        arg_raw.add_argument( "--dna_vcf", dest = "str_dna_vcf", required=True, help = "Input DNA VCF file." )

        # MAF file
        arg_raw.add_argument( "--maf", dest = "str_maf_file", default=None, help = "The maf file used to compare as truth to DNA and RNA vcfs." )
        # MAF file sample id
        arg_raw.add_argument( "--id", dest = "str_maf_sample_id", default=None, help = "The id of the sample in the MAF file." )

        # RNA BAM
        arg_raw.add_argument( "--rna_bam", dest = "str_rna_bam", required=True, help = "Input RNA BAM file." )
        # RNA VCF
        arg_raw.add_argument( "--rna_vcf", dest = "str_rna_vcf", required=True, help = "Input RNA VCF file." )

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """

        * return : List of commands to be ran
                 : list of command objects
        """
        # Check to make sure directories exit
        STR_DEPTH_DIR = os.path.join( args_parsed.str_file_base, self.str_depth_dir )
        STR_LOG_DIR = os.path.join( args_parsed.str_file_base, self.str_log_dir )
        STR_TAB_DIR = os.path.join( args_parsed.str_file_base, self.str_tab_dir )
        STR_FILTERED_VCFS = os.path.join( args_parsed.str_file_base, self.str_filtered_vcf )

        # Get the sample pairing file for the study
        lcmd_commands = []

        # Files for pipeline
        str_RNA_depth_file = os.path.splitext( os.path.basename( args_parsed.str_rna_bam ) )[ 0 ] + "_bam.depth"
        str_RNA_depth_file = os.path.join( STR_DEPTH_DIR, str_RNA_depth_file )
        str_RNA_depth_compressed_file = os.path.splitext( str_RNA_depth_file )[ 0 ] + "_depth.gz"
        str_DNA_depth_file = os.path.splitext( os.path.basename( args_parsed.str_dna_bam ) )[ 0 ] + "_bam.depth"
        str_DNA_depth_file = os.path.join( STR_DEPTH_DIR, str_DNA_depth_file )
        str_DNA_depth_compressed_file = os.path.splitext( str_DNA_depth_file )[ 0 ] + "_depth.gz"
        str_filtered_dna_vcf = os.path.join( STR_FILTERED_VCFS, os.path.basename( os.path.splitext( args_parsed.str_dna_vcf )[ 0 ] + "_snp.vcf" ) )
        str_filtered_rna_vcf = os.path.join( STR_FILTERED_VCFS, os.path.basename( os.path.splitext( args_parsed.str_rna_vcf )[ 0 ] + "_snp.vcf" ) )
        str_file_maf_dna_tab = os.path.join( STR_TAB_DIR, str_sample_key + "_maf_dna.tab" )
        str_log_maf_dna = os.path.join( STR_LOG_DIR, str_sample_key + "_maf_dna.log" )
        str_file_maf_rna_tab = os.path.join( STR_TAB_DIR, str_sample_key + "_maf_rna.tab" )
        str_log_maf_rna = os.path.join( STR_LOG_DIR, str_sample_key + "_maf_rna.log" )
        str_file_dna_rna_tab = os.path.join( STR_TAB_DIR, str_sample_key + "_dna_rna.tab" )
        str_log_dna_rna = os.path.join( STR_LOG_DIR, str_sample_key + "_dna_rna.log" )

        # Create depth file RNA
        lcmd_commands.append( Command.Command( str_cur_command = "samtools depth " + args_parsed.str_rna_bam + " > " + str_RNA_depth_file,
                                               lstr_cur_dependencies = [ args_parsed.str_rna_bam ],
                                               lstr_cur_products = [ str_RNA_depth_file ] ) )
        lcmd_commands.append( Command.Command( str_cur_command = "gzip -c " + str_RNA_depth_file + " > " + str_RNA_depth_compressed_file,
                                               lstr_cur_dependencies = [ str_RNA_depth_file ],
                                               lstr_cur_products = [ str_RNA_depth_compressed_file ] ) )

        # Create and Compress depth file DNA
        # Do not perform command if the file is the same as the RNA
        if not args_parsed.str_dna_bam == args_parsed.str_rna_bam:
            lcmd_commands.append( Command.Command( str_cur_command = "samtools depth " + args_parsed.str_dna_bam + " > " + str_DNA_depth_file,
                                                   lstr_cur_dependencies = [ args_parsed.str_dna_bam ],
                                                   lstr_cur_products = [ str_DNA_depth_file ] ) )
            lcmd_commands.append( Command.Command( str_cur_command = "gzip -c " + str_DNA_depth_file + " > " + str_DNA_depth_compressed_file,
                                                   lstr_cur_dependencies = [ str_DNA_depth_file ],
                                                   lstr_cur_products = [ str_DNA_depth_compressed_file ] ) )

        # Filter RNA and DNA VCF files to SNPS
        # SOFT Filter DNA VCF down to snps
        str_filtered_dna_vcf_command = os.path.join( self.str_src_dir, "reduce_vcf_to_snps.py " ) + args_parsed.str_dna_vcf + " " + str_filtered_dna_vcf
        lcmd_commands.append( Command.Command( str_cur_command = str_filtered_dna_vcf_command,
                                               lstr_cur_dependencies = [ args_parsed.str_dna_vcf ],
                                               lstr_cur_products = [ str_filtered_dna_vcf ] ) )

        # SOFT Filter RNA VCF down to snps
        str_filtered_rna_vcf_command = os.path.join( self.str_src_dir, "reduce_vcf_to_snps.py " ) + args_parsed.str_rna_vcf + " " + str_filtered_rna_vcf
        lcmd_commands.append( Command.Command( str_cur_command = str_filtered_rna_vcf_command,
                                               lstr_cur_dependencies = [ args_parsed.str_rna_vcf ],
                                               lstr_cur_products = [ str_filtered_rna_vcf ] ) )

        # Create tab files
        # Tab files for vcf files that include if common or somatic
        if not args_parsed.str_maf_file is None:
            str_cmd_tab_MAF_DNA = os.path.join( self.str_src_dir, "vcfs_to_snp_calls_tab.py" ) + " --maf_reference " + args_parsed.str_maf_file + " --tumor " + args_parsed.str_maf_sample_id + " --count_reference " + str_DNA_depth_compressed_file + " --vcf " + str_filtered_dna_vcf + " --count " + str_DNA_depth_compressed_file + " " + str_file_maf_dna_tab + " > " + str_log_maf_dna
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_tab_MAF_DNA,
                                                   lstr_cur_dependencies = [ args_parsed.str_maf_file, str_DNA_depth_compressed_file, str_filtered_dna_vcf ],
                                                   lstr_cur_products = [ str_file_maf_dna_tab, str_log_maf_dna ] ) )

            str_cmd_tab_MAF_RNA = os.path.join( self.str_src_dir, "vcfs_to_snp_calls_tab.py" ) + " --maf_reference " + args_parsed.str_maf_file + " --tumor " + args_parsed.str_maf_file + " --count_reference " + str_DNA_depth_compressed_file + " --vcf " + str_filtered_rna_vcf + " --count " + str_RNA_depth_compressed_file + " " + str_file_maf_rna_tab + " > " + str_log_maf_rna
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_tab_MAF_RNA,
                                                   lstr_cur_dependencies = [ args_parsed.str_maf_file, str_DNA_depth_compressed_file, str_RNA_depth_compressed_file, str_filtered_rna_vcf ],
                                                   lstr_cur_products = [ str_file_maf_rna_tab, str_log_maf_rna ] ) )

        str_cmd_tab_DNA_RNA = os.path.join( self.str_src_dir, "vcfs_to_snp_calls_tab.py" ) + " --vcf_reference " + str_filtered_dna_vcf + " --count_reference " + str_DNA_depth_compressed_file + " --vcf " + str_filtered_rna_vcf + " --count " + str_RNA_depth_compressed_file + " " + str_file_dna_rna_tab + " > " + str_log_dna_rna
        lcmd_commands.append( Command.Command( str_cur_command = str_cmd_tab_DNA_RNA,
                                               lstr_cur_dependencies = [ str_filtered_dna_vcf, str_DNA_depth_compressed_file, str_RNA_depth_compressed_file, str_filtered_rna_vcf ],
                                               lstr_cur_products = [ str_file_dna_rna_tab, str_log_dna_rna ] ) )

        return lcmd_commands

if __name__ == "__main__":

    # Run pipeline
    RNASEQ_mutation_make_tabs().func_run_pipeline()
