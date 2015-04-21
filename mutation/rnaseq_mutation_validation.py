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

# Sample study keys
STR_CMDS = "CMDS"
STR_DNA_LEFT = "DNA_LEFT"
STR_DNA_RIGHT = "DNA_RIGHT"
STR_MAF_SAMPLE_NAME = "MAF_SAMPLE"
STR_RNA_LEFT = "RNA_LEFT"
STR_RNA_RIGHT = "RNA_RIGHT"
# STR_SYNTHETIC_LEFT = "SYN_LEFT"
# STR_SYNTHETIC_RIGHT = "SYN_RIGHT"
STR_SYNTHETIC_VCF = "SYN_VCF"

# Paired sample file indices
I_SYNTHETIC_VCF = 5
I_RNA_RIGHT = 4
I_RNA_LEFT = 3
I_DNA_RIGHT = 2
I_DNA_LEFT = 1
I_MAF_SAMPLE = 0

class RNASEQ_mutation_validation( ParentScript.ParentScript ):

    def __init__( self ):

        # File structure
        self.str_bams_dir = "alignment"
        self.str_conf = "conf"
        self.str_figure_dir = "figure"
        self.str_log_dir = "log"
        self.str_sample_files_dir = "sample_files"
        self.str_tab_dir = "tab"
        self.str_truth_runs = "truth_runs"
        self.str_figure_dna = os.path.join( self.str_figure_dir, "dna" )
        self.str_figure_dna_rna = os.path.join( self.str_figure_dir, "dna_rna" )
        self.str_figure_rna = os.path.join( self.str_figure_dir, "rna" )
        self.str_figure_maf = os.path.join( self.str_figure_dir, "maf" )
        self.str_figure_syn = os.path.join( self.str_figure_dir, "syn" )
        self.str_filtered_vcf = "vcf_filtered"

    def parse_paired_sample_file( self, str_paired_file ):
        """
        Parse a file holding the RNA / DNA smaple pairings.

        * return : A dictionary of pairing information stored by RNA sample
                 : Dict { str_rna_sample : { STR_DNA_SAMPLE : str_dna_sample } }
        """

        # Expecting the file to be
        ## mutect_tumor_id  dna_fastq_left  dna_fastq_right rna_fastq_left  rna_fastq_right synthetic_bam synthetic_vcf
        dict_return = {}
        with open( str_paired_file, "r" ) as hndl_paired_file:
          for lstr_line in csv.reader( hndl_paired_file, delimiter = "\t" ):
            if lstr_line[ 0 ][ 0 ] == "#":
              continue
            # Replace NA or None for None to be consistent
            lstr_line = [ None if str_line_token.lower() in ["none", "na" ] else str_line_token for str_line_token in lstr_line ]

            dict_return[ lstr_line[ I_RNA_LEFT ]] = { STR_DNA_RIGHT : lstr_line[ I_DNA_RIGHT ],
                                                   STR_DNA_LEFT : lstr_line[ I_DNA_LEFT ],
                                                   STR_RNA_LEFT : lstr_line[ I_RNA_LEFT ],
                                                   STR_RNA_RIGHT : lstr_line[ I_RNA_RIGHT ],
                                                   STR_MAF_SAMPLE_NAME : lstr_line[ I_MAF_SAMPLE ],
                                                   STR_SYNTHETIC_VCF : lstr_line[ I_SYNTHETIC_VCF ] }
        return dict_return

    def func_update_arguments( self, arg_raw ):
        """
        Add in custom arguments for the mutation validation pipeline.
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "RNASEQ_mutation_validation"
        # Alignment config file. RNASEQ alignment for all variant calling.
        arg_raw.add_argument( "--align_conf", dest = "str_align_run_config", required = True, help = "Alignment config for all variant calling." )

        # Annotation config file for all runs
        arg_raw.add_argument( "--annot_conf", dest = "str_annot_config", required = True, help = "Annotation config file for all samples." )

        # Conf file, one run or all samples per conf file will be ran
        arg_raw.add_argument( "--call_conf", dest = "lstr_run_configs", required = True, action = "append", help = "Config file for the call variants part of the pipeline. Can be used multiple times, each conf file is a setting for a study that is ran on ALL samples and validated as a seperate group." )

        # Number of threads
        arg_raw.add_argument( "--jobs", dest = "i_jobs", default = 1, action = "store", help = "Max number of jobs to run at a time." )

        # Truth calling run conf happends once (on all DNA samples) per validation run.
        arg_raw.add_argument( "--truth_call_conf", dest = "str_truth_run_config", default = None, action = "store", help = "Config file to call variants on the truth data." )

        # File that pairs RNA and DNA seq files
        arg_raw.add_argument( "--sample_file", dest = "str_sample_pairing_file", required = True, help = "A file that contains pairing of RNA and DNA Bam and Vcf files." )

        # MAF file
        # dest = str_maf_file
        arg_raw.add_argument( "--maf", dest = "str_maf_file", default=None, help = "The maf file used to compare DNA dn RNA vcfs." )

        # Reduced dbSNP vcf file
        # dest = str_reduced_dbsnp_vcf
        arg_raw.add_argument( "--ref_vcf", dest = "str_reduced_dbsnp_vcf", required = True, help = "A reference vcf file." )
        arg_raw.add_argument("--ref_gtf", dest = "str_annotation_gtf", required = True, help = "A reference gtf file." )

        # Comma delimited string of key mutations to look at occurence in samples
        # dest = args_parsed.str_key_mutations
        arg_raw.add_argument( "--key", dest = "str_key_mutations", required = True, help = "Comma delimited string of key mutations (gene names) to look at occurence in samples." )

        # Run config file for making tab files
        arg_raw.add_argument( "--tab_conf_dir", dest = "str_tab_config_dir", action = "store", help = "Directory to find run config file for making tab files" )

    def func_convert_fastq_left_rna_bam( self, str_fastq_left, str_output_dir ):
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, self.str_bams_dir, "samples", str_file_base, "COMMON_ALIGNMENT", "star_align_2_reads_"+str_file_base+"_fq_P_qtrim_fq","Aligned.out.bam" )
    def func_convert_fastq_left_rna_vcf( self, str_fastq_left, str_output_dir ):
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, self.str_bams_dir, "samples", str_file_base, "RNASEQ_MUTATION_PREMADE", "variants_annotated.vcf.gz" )
    def func_convert_fastq_left_dna_bam( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "dna.bam" )
    def func_convert_fastq_left_dna_vcf( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "dna_snp.vcf" )
    def func_convert_fastq_maf_dna_tab( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "maf_dna.tab" )
    def func_convert_fastq_maf_rna_tab( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "maf_rna.tab" )
    def func_convert_fastq_dna_rna_tab( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "dna_rna.tab" )
    def func_convert_fastq_syn_rna_tab( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "rna.tab" )
    def func_convert_fastq_rna_depth( self, str_fastq_left, str_output_dir ):
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, self.str_bams_dir, "samples", str_file_base, "COMMON_ALIGNMENT", "Aligned.out.bam.depth.gz" )
    def func_convert_fastq_dna_depth( self, str_fastq_left, str_output_dir ):
        return os.path.join( str_output_dir, "dna.depth" )

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """

        * return : List of commands to be ran
                 : List of command objects
        """

        # Get the sample pairing file for the study
        dict_sample_study = self.parse_paired_sample_file( str_paired_file = args_parsed.str_sample_pairing_file )

        # Commands to run
        lcmd_commands = []
        lcmd_commands_run = []

        # Make needed directory structure
        if not cur_pipeline.func_mkdirs( [  os.path.join( args_parsed.str_file_base, str_file )
                                          for str_file in [ self.str_bams_dir, self.str_figure_dir, self.str_log_dir, 
                                                            self.str_sample_files_dir, self.str_tab_dir, self.str_truth_runs,
                                                            self.str_filtered_vcf, self.str_figure_rna, self.str_figure_dna,
                                                            self.str_figure_dna_rna, self.str_figure_maf ]] ):
            exit( 2 )

        #########################################
        # Prep secondary files if needed.
        # SOFT Filter the DBSNP vcf file to SNP
        #########################################
        str_filtered_dbsnp_vcf = os.path.join( args_parsed.str_file_base, self.str_filtered_vcf, os.path.splitext( os.path.basename( args_parsed.str_reduced_dbsnp_vcf ) )[ 0 ] + "_snp.vcf" )
        str_filtered_dbsnp_vcf_command = os.path.join( "reduce_vcf_to_snps.py --reference " ) + args_parsed.str_reduced_dbsnp_vcf + " " + str_filtered_dbsnp_vcf
        lcmd_commands_run.append( Command.Command( str_cur_command = str_filtered_dbsnp_vcf_command,
                                               lstr_cur_dependencies = [ args_parsed.str_reduced_dbsnp_vcf ],
                                               lstr_cur_products = [ str_filtered_dbsnp_vcf ] ) )

        ###############################
        # Run the truth runs if needed
        ###############################
        lstr_truth_call_products = []
        if args_parsed.str_truth_run_config:
            # Make DNA Seq / truth samples.txt
            str_truth_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "truth_samples.txt" )
            llstr_truth_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_DNA_LEFT ] ) )[0], dict_sample[ STR_DNA_LEFT], dict_sample[ STR_DNA_RIGHT ] ]
              for dict_sample in dict_sample_study.values() if dict_sample[ STR_DNA_LEFT ] and dict_sample[ STR_DNA_RIGHT ] ]
            lstr_truth_samples = [ "\t".join( lstr_files ) for lstr_files in llstr_truth_samples ]
            if len( lstr_truth_samples ) > 0:
                if not args_parsed.f_Test:
                    with open( str_truth_sample_file, "w" ) as hndl_truth_samples:
                        hndl_truth_samples.write( "\n".join( lstr_truth_samples ) )
                # Make dep / product names
                lstr_truth_call_dependencies = [ [ lstr_sample[ 1 ], lstr_sample[2] ] for lstr_sample in llstr_truth_samples ]
                lstr_truth_call_dependencies = [ str_dep for lstr_subgroup in lstr_truth_call_dependencies for str_dep in lstr_subgroup ]
                lstr_truth_call_products = [ [ self.func_convert_fastq_left_dna_bam( lstr_file[ 1 ], args_parsed.str_file_base ), self.func_convert_fastq_left_dna_vcf ( lstr_file[ 1 ], args_parsed.str_file_base ) ] for lstr_file in llstr_truth_samples ]
                lstr_truth_call_products = [ str_prod for lstr_subgroup in lstr_truth_call_products for str_prod in lstr_subgroup ]

                # Run dna samples
                str_cmd_truth_call_variants = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", args_parsed.str_truth_run_config, "--reads_list_file", str_truth_sample_file, 
                                                "--project_base_dir", os.path.join( args_parsed.str_file_base, self.str_truth_runs ),
                                                "--memory 35 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_truth_call_variants,
                                                   lstr_cur_dependencies = lstr_truth_call_dependencies + [ args_parsed.str_annot_config, args_parsed.str_truth_run_config ],
                                                   lstr_cur_products = lstr_truth_call_products ) )

        # Make sample file for rna seq and variant calling from bam
        str_rna_seq_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "rnaseq_samples.txt" )
        str_variant_calling_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "variant_calling_samples.txt" )
        llstr_rna_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0], dict_sample[ STR_RNA_LEFT], dict_sample[ STR_RNA_RIGHT ] ]
                                  for dict_sample in dict_sample_study.values() if dict_sample[ STR_RNA_LEFT ] and dict_sample[ STR_RNA_RIGHT ] ]
        llstr_bam_call_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0], self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT], args_parsed.str_file_base )]
                                  for dict_sample in dict_sample_study.values() if dict_sample[ STR_RNA_LEFT ] ]
        lstr_rna_samples = [ "\t".join( lstr_files ) for lstr_files in llstr_rna_samples ]
        lstr_bam_call_samples = [ "\t".join( lstr_files ) for lstr_files in llstr_bam_call_samples ]
        if len( lstr_rna_samples ) > 0:
            # Make rnaseq alignment sample file
            if not args_parsed.f_Test:
                with open( str_rna_seq_sample_file, "w" ) as hndl_rna_samples:
                    hndl_rna_samples.write( "\n".join( lstr_rna_samples ) )
            # Make variant calling sample file
            if len( lstr_bam_call_samples ) > 0 and not args_parsed.f_Test:
                with open( str_variant_calling_sample_file, "w" ) as hndl_bam_samples:
                    hndl_bam_samples.write( "\n".join( lstr_bam_call_samples ) )
            # Make dep / product names
            lstr_bam_dependencies = [ [ lstr_sample[ 1 ], lstr_sample[2] ] for lstr_sample in llstr_rna_samples ]
            lstr_bam_dependencies = [ str_dep for lstr_subgroup in lstr_bam_dependencies for str_dep in lstr_subgroup ]
            lstr_bam_products = [ self.func_convert_fastq_left_rna_bam( lstr_sample[ 1 ], args_parsed.str_file_base ) for lstr_sample in llstr_rna_samples ]
            lstr_call_dependencies = lstr_bam_products
            lstr_call_products = [ self.func_convert_fastq_left_rna_vcf( lstr_sample[ 1 ], args_parsed.str_file_base ) for lstr_sample in llstr_rna_samples ]

            ########################
            # Align RNA Seq bams
            ########################
            str_cmd_make_bam = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                            "--run_conf", args_parsed.str_align_run_config, "--reads_list_file", str_rna_seq_sample_file,
                                            "--project_base_dir", os.path.join( args_parsed.str_file_base, self.str_bams_dir),
                                            "--memory 35 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
            lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_make_bam,
                                                   lstr_cur_dependencies = lstr_bam_dependencies + [ args_parsed.str_annot_config, args_parsed.str_align_run_config ],
                                                   lstr_cur_products = lstr_bam_products ) )
           
            # For each conf file
            for str_call_run_conf in args_parsed.lstr_run_configs:
                # Current project directory
                str_current_project_dir = os.path.join( args_parsed.str_file_base, os.path.splitext( os.path.basename( str_call_run_conf ) )[ 0 ] )
                if not cur_pipeline.func_mkdirs( [ str_current_project_dir ] ):
                    exit( 3 )
                # Track files made
                lstr_generated_dna_vcf_snps = []
                lstr_generated_rna_vcf_snps = []
                lstr_generated_maf_rna_tab = []
                lstr_generated_maf_dna_tab = []
                lstr_generated_dna_rna_tab = []
                lstr_generated_syn_rna_tab = []

                #############################
                # Call variants off of a bam
                #############################
                str_cmd_call_variants = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", str_call_run_conf, "--reads_list_file", str_variant_calling_sample_file, 
                                                "--project_base_dir", str_current_project_dir, "--memory 35 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                lstr_call_products = [ os.path.join( str_current_project_dir, "samples", os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0],
                                       "RNASEQ_MUTATION_PREMADE", os.path.basename( str_file ) ) for str_file in lstr_call_products ]
                lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_call_variants,
                                                   lstr_cur_dependencies = lstr_call_dependencies + [ args_parsed.str_annot_config, str_call_run_conf ],
                                                   lstr_cur_products = lstr_call_products ) )

                # Make tabs for MAF and RNASEQ Data and/or DNASEQ and RNASEQ Data
                if args_parsed.str_truth_run_config or args_parsed.str_maf_file:
                    # Make the sample text file for tabs
                    # Content will be decided depending on the type of run
                    str_tab_sample_file = os.path.join( str_current_project_dir, "tab_samples.txt" )
                    # Make products and dependencies
                    # Content will be decided depending on the type of run
                    lstr_tab_dependencies = []
                    lstr_tab_products = []
                    llstr_tab_samples = []
                    lstr_link_dependencies = []
                    lstr_link_products = []
 
                    # Default just MAF conf, will be written over if another conf is needed.
                    str_tab_run_conf_file = None
                    lstr_symbolic_link_command = []
                    
                    if args_parsed.str_truth_run_config:
                        if args_parsed.str_maf_file:
                            str_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, "run_maf_dnaseq_tabs.conf" )
                            lstr_tab_dependencies.append( args_parsed.str_maf_file )
                            for dict_sample in dict_sample_study.values():
                                if not dict_sample[ STR_DNA_LEFT ] or not dict_sample[ STR_DNA_RIGHT ]:
                                    continue
                                lstr_tab_products.append( self.func_convert_fastq_maf_dna_tab( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_products.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_products.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_dna_bam( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
#                                lstr_generated_dna_vcf_snps.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
#                                lstr_generated_rna_vcf_snps.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_generated_maf_dna_tab.append( self.func_convert_fastq_maf_dna_tab( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_generated_maf_rna_tab.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_generated_dna_rna_tab.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )

                                str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0]
                                lstr_symbolic_link_command.append( "ln -sf " + args_parsed.str_maf_file + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".maf" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + "_reference.vcf.gz" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_left_rna_snp_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".vcf.gz" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_dna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + "_reference.depth.gz" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_rna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".depth.gz" ) )
                                lstr_link_products = [ os.path.join( str_current_project_dir, str_cur_tab_sample + str_ext ) for str_ext in ["_reference.vcf",".vcf.gz","_reference.depth.gz",".depth.gz",".maf"] ]

                        else:
                            str_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, "run_dnaseq_tabs.conf" )
                            for dict_sample in dict_sample_study.values():
                                if not dict_sample[ STR_DNA_LEFT ] or not dict_sample[ STR_DNA_RIGHT ]:
                                    continue
                                lstr_tab_products.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_dna_bam( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_tab_dependencies.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
#                                lstr_generated_dna_vcf_snps.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
#                                lstr_generated_rna_vcf_snps.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                                lstr_generated_dna_rna_tab.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )

                                str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0]
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + "_reference.vcf.gz" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".vcf.gz" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_dna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + "_reference.depth.gz" ) )
                                lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_rna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".depth.gz" ) )
                                lstr_link_products = [ os.path.join( str_current_project_dir, str_cur_tab_sample + str_ext ) for str_ext in ["_reference.vcf",".vcf.gz","_reference.depth.gz",".depth.gz"] ]
                    elif args_parsed.str_maf_file:
                        str_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, "run_maf_tabs.conf" )
                        for dict_sample in dict_sample_study.values():
                            if not dict_sample[ STR_MAF_SAMPLE_NAME ]:
                                continue
                            lstr_tab_products.append( self.func_convert_fastq_maf_dna_tab( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                            lstr_tab_products.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_left_dna_bam( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
#                            lstr_generated_dna_vcf_snps.append( self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
#                            lstr_generated_rna_vcf_snps.append( self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_generated_maf_rna_tab.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )

                            str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0]
                            lstr_symbolic_link_command.append( "ln -sf " + args_parsed.str_maf_file + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".maf" ) )
                            lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".vcf.gz" ) )
                            lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_rna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".depth.gz" ) )
                            lstr_link_products = [ os.path.join( str_current_project_dir, str_cur_tab_sample + str_ext ) for str_ext in [".maf",".vcf.gz",".depth.gz"] ]

                    # Make the samples file for the tabs
                    for dict_sample in dict_sample_study.values():
                        if not dict_sample[ STR_RNA_LEFT ] or not dict_sample[ STR_RNA_RIGHT ]:
                            continue
                        str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0]
                        llstr_tab_samples.append( "\t".join( [ str_cur_tab_sample,
                                                        os.path.join( str_current_project_dir, str_cur_tab_sample ),
                                                        str( dict_sample[ STR_MAF_SAMPLE_NAME ] ) + "\n" ] ) )
                    if len( llstr_tab_samples ) > 0 and len( lstr_symbolic_link_command ) > 0:
                        if not args_parsed.f_Test:
                            with open( str_tab_sample_file, "w" ) as hndl_tab_samples:
                                hndl_tab_samples.write( "\n".join( llstr_tab_samples ) )

                        # Make symbolic links
                        str_command_symbolic_link = ";".join( lstr_symbolic_link_command )
                        lcmd_commands_run.append( Command.Command( str_cur_command = str_command_symbolic_link,
                                                       lstr_cur_dependencies = [ args_parsed.str_annot_config, str_tab_run_conf_file ],#TODO Add ln dependencies
                                                       lstr_cur_products = lstr_link_products ) )
                        # Make tabs for biological truth data and RNASEQ Data
                        str_cmd_make_tabs = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", str_tab_run_conf_file, "--reads_list_file", str_tab_sample_file,
                                                "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parrallel_procs", str( args_parsed.i_jobs ) ] )
                        lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_make_tabs,
                                                       lstr_cur_dependencies = lstr_tab_dependencies + [ args_parsed.str_annot_config, str_tab_run_conf_file ],
                                                       lstr_cur_products = lstr_tab_products ) )

                # Make synthetic comparisons tab files
                # Make tabs for synthetic truth data and RNASEQ Data
                # Formats: Synthetic
                # ID from RNA left \t Synthetic bam \t Synthetic vcf \t RNA bam \t RNA vcf
                str_synthetic_tab_sample_file = os.path.join( str_current_project_dir, "tab_syn_samples.txt" )
                str_synthetic_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, "run_synthetic_tabs.conf" )
                lstr_syn_tab_dependencies = []
                lstr_syn_tab_products = []
                llstr_synthetic_samples = []
                lstr_symbolic_link_command = []
                for dict_sample in dict_sample_study.values():
                    if not dict_sample[ STR_SYNTHETIC_VCF ]:
                        continue
                    str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0]
                    llstr_synthetic_samples.append( "\t".join( [ str_cur_tab_sample,
                                               os.path.join( str_current_project_dir, str_cur_tab_sample ),
                                               str( dict_sample[ STR_MAF_SAMPLE_NAME ] ) + "\n" ] ) )
                    lstr_syn_tab_dependencies.extend( [ dict_sample[ STR_SYNTHETIC_VCF ],
                                             self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ),
                                             self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) ] )
                    lstr_syn_tab_products.append( self.func_convert_fastq_syn_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                    lstr_generated_syn_rna_tab.append(  self.func_convert_fastq_syn_rna_tab( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) )
                    lstr_symbolic_link_command.append( "ln -sf " + dict_sample[ STR_SYNTHETIC_VCF ] + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + "_reference.vcf" ) )
                    lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".vcf.gz" ) )
                    lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_rna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + "_reference.depth.gz" ) )
                    lstr_symbolic_link_command.append( "ln -sf " + self.func_convert_fastq_rna_depth( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ) + " " + os.path.join( str_current_project_dir, str_cur_tab_sample + ".depth.gz" ) )
                    lstr_link_products = [ os.path.join( str_current_project_dir, str_cur_tab_sample + str_ext ) for str_ext in ["_reference.vcf",".vcf.gz","_reference.depth.gz",".depth.gz"] ]
                if len( llstr_synthetic_samples ) > 0 and len( lstr_symbolic_link_command ) > 0:
                    if not args_parsed.f_Test:
                        with open( str_synthetic_tab_sample_file, "w" ) as hndl_syn_samples:
                            hndl_syn_samples.write( "\n".join( llstr_synthetic_samples ) )
                    str_cmd_make_synthetic_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", str_synthetic_tab_run_conf_file, "--reads_list_file", str_synthetic_tab_sample_file,
                                                "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parrallel_procs", str( args_parsed.i_jobs ) ] )
                    # Make symbolic links
                    str_command_symbolic_link = ";".join( lstr_symbolic_link_command )
                    lcmd_commands_run.append( Command.Command( str_cur_command = str_command_symbolic_link,
                                                               lstr_cur_dependencies = [ args_parsed.str_annot_config, str_tab_run_conf_file ],
                                                               lstr_cur_products = lstr_link_products ) )
                    lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_make_synthetic_tab,
                                                       lstr_cur_dependencies = lstr_syn_tab_dependencies + [ args_parsed.str_annot_config, str_synthetic_tab_run_conf_file ],
                                                       lstr_cur_products = lstr_syn_tab_products ) )

                # Make figures for validation with a study
                dict_validation_commands = self.func_validation_figure_commands( args_parsed = args_parsed, str_cur_project_dir = str_current_project_dir,
                                                                        f_synthetic = not llstr_synthetic_samples is None,
                                                                        lstr_dna_vcfs_snps = lstr_generated_dna_vcf_snps,
                                                                        lstr_rna_vcfs_snps = lstr_generated_dna_vcf_snps,
                                                                        lstr_maf_rna_tab = lstr_generated_maf_rna_tab,
                                                                        lstr_maf_dna_tab = lstr_generated_maf_dna_tab,
                                                                        lstr_dna_rna_tab = lstr_generated_dna_rna_tab,
                                                                        lstr_syn_rna_tab = lstr_generated_syn_rna_tab ) 
                lcmd_commands.extend( dict_validation_commands[ STR_CMDS ] )

        # TODO Add in figures for validation between studies    
        # This includes ROCs of methods side by side so ROCs for N! number of run confs.
        return( lcmd_commands_run )

    def func_validation_figure_commands( self, args_parsed, str_cur_project_dir, f_synthetic, lstr_dna_vcfs_snps = [], lstr_rna_vcfs_snps = [], lstr_maf_rna_tab = [], lstr_maf_dna_tab = [], lstr_dna_rna_tab = [], lstr_syn_rna_tab = [] ):
        """

        * return : List of commands to be ran
                 : list of command objects
        """

        # Collects commands
        lcmd_commands = []

        # Distance metric functions
        str_jaccard_distance = os.path.join( "jaccard_distance.R" )
        str_jaccard_distance_restricted = os.path.join( "jaccard_distance_restricted.R" )

        # Check to make sure directories exit
        STR_FIGURE_DIR = os.path.join( str_cur_project_dir, self.str_figure_dir )
        STR_DNA_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_dna )
        STR_RNA_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_rna )
        STR_MAF_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_maf )
        STR_SYN_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_syn )

        # Plot the rate of mutation over a window.
        # Mutation rates are often local.
        # Would be interesting to filter on # of mutations in a gene given the mutation rates locally.
        if( len( lstr_dna_vcfs_snps ) > 1 and len( lstr_rna_vcfs_snps ) > 1 ):
            # Create genotype matrix and list
            ## DNA and RNA
            str_genotype_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype_matrix.txt" )
            str_cmd_create_genotype_matrix = "vcfs_to_genotype_matrix.py" + " --matrix " + str_genotype_matrix + " " + " ".join( lstr_rna_vcfs_snps + lstr_dna_vcfs_snps )
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_create_genotype_matrix,
                                               lstr_cur_dependencies = lstr_rna_vcfs_snps + lstr_dna_vcfs_snps,
                                               lstr_cur_products = [ str_genotype_matrix ] ) )
            # Visualize genotype matrices
            ## DNA and RNA, restricted by type
            str_genotype_matrix_restricted_pdf = os.path.join( STR_FIGURE_DIR, "sample_dna_rna_genotype_matrix_restricted.pdf" )
            str_genotype_distance_restricted_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype_restricted.dist" )
            str_cmd_genotype_matrix_dna_rna_restricted = " ".join( [ "make_dendrogram_generic.R", "--input_matrix", str_genotype_matrix, 
                                                      "--distance_function", str_jaccard_distance_restricted,
                                                      "--output_pdf", str_genotype_matrix_restricted_pdf,
                                                      "--output_distance_matrix", str_genotype_distance_restricted_matrix ] ) 
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_dna_rna_restricted,
                                               lstr_cur_dependencies = [ str_genotype_matrix, str_jaccard_distance_restricted ],
                                               lstr_cur_products = [ str_genotype_matrix_restricted_pdf, str_genotype_distance_restricted_matrix ] ) )
            ## DNA and RNA not restricted by type
            str_genotype_matrix_pdf = os.path.join( STR_FIGURE_DIR, "sample_dna_rna_genotype_matrix.pdf" )
            str_genotype_distance_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype.dist" )
            str_cmd_genotype_matrix_dna_rna = " ".join( [ "make_dendrogram_generic.R", "--input_matrix", str_genotype_matrix, 
                                                      "--distance_function", str_jaccard_distance,
                                                      "--output_pdf", str_genotype_matrix_pdf,
                                                      "--output_distance_matrix", str_genotype_distance_matrix ] ) 
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_dna_rna,
                                               lstr_cur_dependencies = [ str_genotype_matrix, str_jaccard_distance ],
                                               lstr_cur_products = [ str_genotype_matrix_pdf, str_genotype_distance_matrix ] ) )

        if( len( lstr_rna_vcfs_snps ) > 1 ):
            ## RNA
            str_genotype_rna_matrix = os.path.join( args_parsed.str_file_base, "sample_rna_genotype_matrix.txt" )
            str_cmd_create_genotype_rna_matrix = "vcfs_to_genotype_matrix.py" + " --matrix " + str_genotype_rna_matrix + " " + " ".join( lstr_rna_vcfs_snps )
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_create_genotype_rna_matrix,
                                               lstr_cur_dependencies = lstr_rna_vcfs_snps,
                                               lstr_cur_products = [ str_genotype_rna_matrix ] ) )
        
            # Visualize genotype matrices
            ## RNA
            str_genotype_rna_matrix_pdf = os.path.join( STR_FIGURE_DIR, "sample_rna_genotype_matrix.pdf" )
            str_genotype_rna_distance_matrix = os.path.join( args_parsed.str_file_base, "sample_rna_genotype.dist" )
            str_cmd_genotype_matrix_rna = " ".join( [ os.path.join( "make_dendrogram_generic.R" ), "--input_matrix", str_genotype_rna_matrix, 
                                                      "--distance_function", str_jaccard_distance,
                                                      "--output_pdf", str_genotype_rna_matrix_pdf,
                                                      "--output_distance_matrix", str_genotype_rna_distance_matrix ] )
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_rna,
                                               lstr_cur_dependencies = [ str_genotype_rna_matrix, str_jaccard_distance ],
                                               lstr_cur_products = [ str_genotype_rna_matrix_pdf, str_genotype_rna_distance_matrix ] ) )

        # Explore false positive and negative rates
        # DNA_RNA
        if len( lstr_dna_rna_tab ):
            str_DNA_RNA_figures = os.path.join( STR_FIGURE_DIR, "dna_rna" )
            lstr_DNA_RNA_roc_1 = [ os.path.join( str_DNA_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_dna_rna_tab ]
            lstr_DNA_RNA_roc_2 = [ os.path.join( str_DNA_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_dna_rna_tab ]
            str_compare_DNA_RNA_cmd = "visualize_mutation_depth_tab_files.R" + " -o " + str_DNA_RNA_figures + " -k DNA_RNA --serial_plots " + " ".join( lstr_dna_rna_tab )
            lcmd_commands.append( Command.Command( str_cur_command = str_compare_DNA_RNA_cmd,
                                               lstr_cur_dependencies = lstr_dna_rna_tab,
                                               lstr_cur_products = lstr_DNA_RNA_roc_1 + lstr_DNA_RNA_roc_2 ) )

        # SYNTHETIC
        if len( lstr_syn_rna_tab ):
            str_SYN_RNA_figures = os.path.join( STR_FIGURE_DIR, "syn_rna" )
            lstr_SYN_RNA_roc_1 = [ os.path.join( str_SYN_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_syn_rna_tab ]
            lstr_SYN_RNA_roc_2 = [ os.path.join( str_SYN_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_syn_rna_tab ]
            str_compare_SYN_RNA_cmd = "visualize_mutation_depth_tab_files.R" + " -o " + str_SYN_RNA_figures + " -k SYN_RNA --serial_plots " + " ".join( lstr_syn_rna_tab )
            lcmd_commands.append( Command.Command( str_cur_command = str_compare_SYN_RNA_cmd,
                                               lstr_cur_dependencies = lstr_syn_rna_tab,
                                               lstr_cur_products = lstr_SYN_RNA_roc_1 + lstr_SYN_RNA_roc_2 ) )

        if not args_parsed.str_maf_file is None:
            # MAF_DNA
            str_MAF_DNA_figures = os.path.join( STR_FIGURE_DIR, "maf_dna" )
            lstr_MAF_DNA_roc_1 = [ os.path.join( str_MAF_DNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_maf_dna_tab ]
            lstr_MAF_DNA_roc_2 = [ os.path.join( str_MAF_DNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_maf_dna_tab ]
            str_compare_MAF_DNA_cmd = "visualize_mutation_depth_tab_files.R" + " -o " + str_MAF_DNA_figures + " -k MAF_DNA --serial_plots " + " ".join( lstr_maf_dna_tab )
            lcmd_commands.append( Command.Command( str_cur_command = str_compare_MAF_DNA_cmd,
                                                   lstr_cur_dependencies = lstr_maf_dna_tab,
                                                   lstr_cur_products = lstr_MAF_DNA_roc_1 + lstr_MAF_DNA_roc_2 ) )

            # MAF_RNA
            str_MAF_RNA_figures = os.path.join( STR_FIGURE_DIR, "maf_rna" )
            lstr_MAF_RNA_roc_1 = [ os.path.join( str_MAF_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_maf_rna_tab ]
            lstr_MAF_RNA_roc_2 = [ os.path.join( str_MAF_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_maf_rna_tab ]
            str_compare_MAF_RNA_cmd = "visualize_mutation_depth_tab_files.R" + " -o " + str_MAF_RNA_figures + " -k MAF_RNA --serial_plots " + " ".join( lstr_maf_rna_tab )
            lcmd_commands.append( Command.Command( str_cur_command = str_compare_MAF_RNA_cmd,
                                                   lstr_cur_dependencies = lstr_maf_rna_tab,
                                                   lstr_cur_products = lstr_MAF_RNA_roc_1 + lstr_MAF_RNA_roc_2 ) )

            if not f_synthetic and args_parsed.str_key_mutations:
                # MAF
                str_key_mutations_tab = "".join( [ " -t " + str_tab for str_tab in lstr_maf_dna_tab ] )
                str_key_mutation_MAF_DNA_output_pdf = os.path.join( STR_MAF_FIGURE, "key_mutation_counts.pdf" )
                str_key_mutations_cmd = "tab_to_percent_mutations.py" + " --gtf " + args_parsed.str_annotation_gtf + " --key " + args_parsed.str_key_mutations + " -o " + str_key_mutation_MAF_DNA_output_pdf + str_key_mutations_tab
                lcmd_commands.append( Command.Command( str_cur_command = str_key_mutations_cmd,
                                               lstr_cur_dependencies = lstr_maf_dna_tab + [ args_parsed.str_annotation_gtf ],
                                               lstr_cur_products = [ str_key_mutation_MAF_DNA_output_pdf ] ) )
                # RNA
                str_key_mutations_tab = "".join( [ " -t " + str_tab for str_tab in lstr_maf_rna_tab ] )
                str_key_mutation_MAF_RNA_output_pdf = os.path.join( STR_RNA_FIGURE, "key_mutation_counts.pdf" )
                str_key_mutations_cmd = "tab_to_percent_mutations.py" + " --second --gtf " + args_parsed.str_annotation_gtf + " --key " + args_parsed.str_key_mutations + " -o " + str_key_mutation_MAF_RNA_output_pdf + str_key_mutations_tab
                lcmd_commands.append( Command.Command( str_cur_command = str_key_mutations_cmd,
                                                   lstr_cur_dependencies = lstr_maf_rna_tab + [ args_parsed.str_annotation_gtf ],
                                                   lstr_cur_products = [ str_key_mutation_MAF_RNA_output_pdf ] ) )

        if not f_synthetic and args_parsed.str_key_mutations:
            # Count mutations and plot
            # DNA
            str_key_mutations_tab = "".join( [ " -t " + str_tab for str_tab in lstr_dna_rna_tab ] )
            str_key_mutation_DNA_RNA_output_pdf = os.path.join( STR_DNA_FIGURE, "key_mutation_counts.pdf" )
            str_key_mutations_cmd = "tab_to_percent_mutations.py" + " --gtf " + args_parsed.str_annotation_gtf + " --key " + args_parsed.str_key_mutations + " -o " + str_key_mutation_DNA_RNA_output_pdf + str_key_mutations_tab
            lcmd_commands.append( Command.Command( str_cur_command = str_key_mutations_cmd,
                                               lstr_cur_dependencies = lstr_dna_rna_tab + [ args_parsed.str_annotation_gtf ],
                                               lstr_cur_products = [ str_key_mutation_DNA_RNA_output_pdf ] ) )

        return { STR_CMDS : lcmd_commands }

if __name__ == "__main__":

    # Run pipeline
    RNASEQ_mutation_validation().func_run_pipeline()
