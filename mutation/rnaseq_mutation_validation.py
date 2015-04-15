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
STR_DNA_LEFT = "DNA_LEFT"
STR_DNA_RIGHT = "DNA_RIGHT"
STR_MAF_SAMPLE_NAME = "MAF_SAMPLE"
STR_RNA_LEFT = "RNA_LEFT"
STR_RNA_RIGHT = "RNA_RIGHT"
STR_SYNTHETIC_BAM = "SYN_BAM"
STR_SYNTHETIC_VCF = "SYN_VCF"

# Paired sample file indices
I_SYNTHETIC_VCF = 6
I_SYNTHETIC_BAM = 5
I_RNA_RIGHT = 4
I_RNA_LEFT = 3
I_DNA_RIGHT = 2
I_DNA_LEFT = 1
I_MAF_SAMPLE = 0

class RNASEQ_mutation_validation( ParentScript.ParentScript ):

    def __init__( self ):

        # File structure
        self.str_depth_dir = "depth"
        self.str_figure_dir = "figure"
        self.str_igv = "igv"
        self.str_log_dir = "log"
        self.str_sample_files_dir = "sample_files"
        self.str_src_dir = "src"
        self.str_tab_dir = "tab"
        self.str_truth_runs = "truth_runs"
        self.str_figure_dna = "dna"
        self.str_figure_rna = "rna"
        self.str_figure_maf = "maf"
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
        f_skip_header = True
        with open( str_paired_file, "r" ) as hndl_paired_file:
          for lstr_line in csv.reader( hndl_paired_file, delimiter = "\t" ):
            if f_skip_header:
              f_skip_header = False
              continue
            # Replace NA or None for None to be consistent
            lstr_line = [ None if str_line_token.lower() in ["none", "na" ] else str_line_token for str_line_token in lstr_line ]

            dict_return[ lstr_line[ I_DNA_LEFT ]] = { STR_DNA_RIGHT : lstr_line[ I_DNA_RIGHT ],
                                                   STR_RNA_LEFT : lstr_line[ I_RNA_LEFT ],
                                                   STR_RNA_RIGHT : lstr_line[ I_RNA_RIGHT ],
                                                   STR_MAF_SAMPLE_NAME : lstr_line[ I_MAF_SAMPLE ],
                                                   STR_SYNTHETIC_BAM : lstr_line[ I_SYNTHETIC_BAM ],
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
        # Annotation config file for all runs
        arg_raw.add_argument( "--annot_conf", dest = "str_annot_config", required = True, help = "Annotation config file for all samples." )

        # Conf file, one run or all samples per conf file will be ran
        arg_raw.add_argument( "--call_conf", dest = "lstr_run_configs", required = True, action = "append", nargs = "+", help = "Config file for the call variants part of the pipeline. Can be used multiple times, each conf file is a setting for a study that is ran on ALL samples and validated as a seperate group." )

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

        # Comma delimited string of key mutations to look at occurence in samples
        # dest = args_parsed.str_key_mutations
        arg_raw.add_argument( "--key", dest = "str_key_mutations", required = True, help = "Comma delimited string of key mutations (gene names) to look at occurence in samples." )

        # Run config file for making tab files
        arg_raw.add_argument( "--tab_conf", dest = "str_tab_config", action = "store", help = "Run config file for making tab files" )

    def func_convert_fastq_left_dna_bam( self, str_fastq_left ):
        return str_fastq_left + "_dna.bam"
    def func_convert_fastq_left_dna_vcf( self, str_fastq_left ):
        return str_fastq_left + "_dna.vcf"
    def func_convert_fastq_left_rna_bam( self, str_fastq_left ):
        return str_fastq_left + "_rna.bam"
    def func_convert_fastq_left_rna_vcf( self, str_fastq_left ):
        return str_fastq_left + "_rna.vcf"
    def func_convert_fastq_maf_dna_tab( self, str_fastq_left ):
        return str_fastq_left + "_maf_dna.tab"
    def func_convert_fastq_maf_rna_tab( self, str_fastq_left ):
        return str_fastq_left + "_maf_rna.tab"
    def func_convert_fastq_dna_rna_tab( self, str_fastq_left ):
        return str_fastq_left +  "_dna_rna.tab"
    def func_convert_fastq_syn_rna_tab( self, str_fastq_left ):
        return str_fastq_left + "_syn_rna.tab"
    def func_convert_fastq_rna_snp_vcf( self, str_fastq_left ):
        return str_fastq_left + "_rna_snp.tab"
    def func_convert_fastq_dna_snp_vcf( self, str_fastq_left ):
        return str_fastq_left + "_dna_snp.vcf"

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """

        * return : List of commands to be ran
                 : List of command objects
        """

        # Get the sample pairing file for the study
        dict_sample_study = self.parse_paired_sample_file( str_paired_file = args_parsed.str_sample_pairing_file )

        # Commands to run
        lcmd_commands = []

        # Prep secondary files if needed.
        # SOFT Filter the DBSNP vcf file to SNP
        str_filtered_dbsnp_vcf = os.path.join( STR_FILTERED_VCFS, os.path.splitext( os.path.basename( args_parsed.str_reduced_dbsnp_vcf ) )[ 0 ] + "_snp.vcf" )
        str_filtered_dbsnp_vcf_command = os.path.join( self.str_src_dir, "reduce_vcf_to_snps.py --reference " ) + args_parsed.str_reduced_dbsnp_vcf + " " + str_filtered_dbsnp_vcf
        lcmd_commands.append( Command.Command( str_cur_command = str_filtered_dbsnp_vcf_command,
                                               lstr_cur_dependencies = [ args_parsed.str_reduced_dbsnp_vcf ],
                                               lstr_cur_products = [ str_filtered_dbsnp_vcf ] ) )

        # Run the truth runs if needed
        if args_parsed.str_truth_run_config:
            # Make DNA Seq / truth samples.txt
            str_truth_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "truth_samples.txt" )
            llstr_truth_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_DNA_LEFT ] ) )[0], dict_sample[ STR_DNA_LEFT], dict_sample[ STR_DNA_RIGHT ] + "\n" ]
              for dict_sample in dict_sample_study.items() if dict_sample[ STR_DNA_LEFT ] and dict_sample[ STR_DNA_RIGHT ] ]
            with open( str_truth_sample_file, "w" ) as hndl_truth_samples:
                hndl_truth_samples.write( llstr_truth_samples )
            # Make dep / product names
            lstr_truth_call_dependencies = [ [ lstr_sample[ 1 ], lstr_sample[2][:-1] ] for lstr_sample in llstr_truth_samples ]
            lstr_truth_call_dependencies = [ str_dep for str_dep in lstr_subgroup for lstr_subgroup in lstr_truth_call_dependencies ]
            lstr_truth_call_products = [ [ self.func_convert_fast_left_dna_vcf( lstr_sample[ 1 ] ), self.func_convert_fast_left_dna_bam( lstr_sample[ 1 ] ) ] for lstr_sample  in llstr_truth_samples ]
            lstr_truth_call_products = [ str_prod for lstr_subgroup in lstr_truth_call_dependencies for str_prod in lstr_subgroup ]

            # Run dna samples
            str_cmd_truth_call_variants = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", args_parsed.str_truth_run_config, "--reads_list_file", str_truth_sample_file, 
                                                "--project_base_dir", os.path.join( args_parsed.str_file_base, self.str_truth_runs ),
                                                "--memory 35 --run_on_grid --num_parallel_procs", args_parsed.i_jobs ] )
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_truth_call_variants,
                                                   lstr_cur_dependencies = lstr_truth_call_dependencies,
                                                   lstr_products = lstr_truth_call_products ) )

        # Make sample file for rna seq
        str_rna_seq_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "rnaseq_samples.txt" )
        llstr_rna_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0], dict_sample[ STR_RNA_LEFT], dict_sample[ STR_RNA_RIGHT ] + "\n" ]
                                  for dict_sample in dict_sample_study.items() if dict_sample[ STR_RNA_LEFT ] and dict_sample[ STR_RNA_RIGHT ] ]
        with open( str_rna_seq_sample_file, "w" ) as hndl_rna_samples:
            hndl_rna_samples.write( llstr_rna_samples )
        # Make dep / product names
        lstr_call_dependencies = [ [ lstr_sample[ 1 ], lstr_sample[2][:-1] ] for lstr_sample in llstr_rna_samples ]
        lstr_call_dependencies = [ str_dep for lstr_subgroup in lstr_call_dependencies for str_dep in lstr_subgroup ]
        lstr_call_products = [ [ self.func_convert_fast_left_rna_vcf( lstr_sample[ 1 ] ), self.func_convert_fast_left_rna_bam( lstr_sample[ 1 ] ) ] for lstr_sample  in llstr_rna_samples ]
        lstr_call_products = [ str_prod for lstr_subgroup in lstr_call_dependencies for str_prod in lstr_subgroup ]

        # For each conf file
        for str_call_run_conf in lstr_run_configs:
            # Current project directory
            str_current_project_dir = os.path.join( arg_parsed.str_file_base, os.path.splitext( os.path.basename( str_call_run_conf ) )[ 0 ] )

            # Track files made
            lstr_generated_dna_vcf_snps = []
            lstr_generated_rna_vcf_snps = []
            lstr_generated_maf_rna_tab = []
            lstr_generated_maf_rna_tab = []
            lstr_generated_dna_rna_tab = []
            lstr_generated_syn_rna_tab = []

            # Call variantss
            str_cmd_call_variants = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annotation_config,
                                                "--run_conf", str_call_run_conf, "--reads_list_file", str_rna_seq_sample_file, 
                                                "--project_base_dir", str_current_project_dir, "--memory 35 --run_on_grid --num_parallel_procs", args_parsed.i_jobs ] )
            lstr_call_dependencies = lstr_call_dependencies 
            lstr_call_products = [ os.path.join( str_cur_project_dir, os.path.basename( str_file ) ) for str_file in lstr_call_products ]
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_call_variants,
                                                   lstr_cur_dependencies = lstr_call_dependencies,
                                                   lstr_products = lstr_call_products ) )

            # Make tabs for MAF and RNASEQ Data and/or DNASEQ and RNASEQ Data
            if args_parsed.str_truth_run_config or args_parsed.str_maf_file:
                # Make the sample text file for tabs
                # Content will be decided depending on the type of run
                str_tab_sample_file = os.path.join( str_current_project_dir, "tab_samples.txt" )
                # Make products and dependencies
                # Content will be decided depending on the type of run
                lstr_tab_dependencies = []
                lstr_tab_products = []
 
                # Default just MAF conf, will be written over if another conf is needed.
                str_tab_run_conf_file = "run_maf_tabs.conf"
                if args_parsed.str_truth_run_config:
                    if args_parsed.str_maf_file:
                        str_tab_run_conf_file = "run_maf_dnaseq_tabs.conf"
                        lstr_tab_dependencies.append( args_parsed.str_maf_file )
                        for dict_sample in dict_sample_study.items():
                            if dict_sample[ STR_DNA_LEFT ] and dict_sample[ STR_DNA_RIGHT ]:
                                continue
                            lstr_tab_products.append( self.func_convert_fastq_maf_dna_tab( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_tab_products.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_tab_products.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_dna_bam( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_rna_bam( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_generated_dna_vcf_snps.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_generated_rna_vcf_snps.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_generated_maf_dna_tab.append( self.func_convert_fastq_maf_dna_tab( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_generated_maf_rna_tab.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_generated_dna_rna_tab.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                            llstr_tab_samples.append( [ dict_sample[ STR_MAF_SAMPLE_NAME ],
                                                        self.func_convert_fastq_dna_bam( dict_sample[ STR_DNA_LEFT ] ),
                                                        self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ),
                                                        self.func_convert_fastq_rna_bam( dict_sample[ STR_RNA_LEFT ] ),
                                                        self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) + "\n" ] )

                    else:
                        str_tab_run_conf_file = "run_dnaseq_tabs.conf"
                        for dict_sample in dict_sample_study.items():
                            if dict_sample[ STR_DNA_LEFT ] and dict_sample[ STR_DNA_RIGHT ]:
                                continue
                            lstr_tab_products.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_dna_bam( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_rna_bam( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_tab_dependencies.append( self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_generated_dna_vcf_snps.append( self.func_convert_fastq_left_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
                            lstr_generated_rna_vcf_snps.append( self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                            lstr_generated_dna_rna_tab.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                            llstr_tab_samples.append( [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0],
                                                        self.func_convert_fastq_dna_bam( dict_sample[ STR_DNA_LEFT ] ),
                                                        self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ),
                                                        self.func_convert_fastq_rna_bam( dict_sample[ STR_RNA_LEFT ] ),
                                                        self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) + "\n" ] )

                else:
                    for dict_sample in dict_sample_study.items():
                        if dict_sample[ STR_RNA_LEFT ] and dict_sample[ STR_RNA_RIGHT ]:
                            continue
                        lstr_tab_products.append( self.func_convert_fastq_maf_dna_tab( dict_sample[ STR_DNA_LEFT ] ) )
                        lstr_tab_products.append( self.func_convert_fastq_maf_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                        lstr_tab_dependencies.append( self.func_convert_fastq_dna_bam( dict_sample[ STR_DNA_LEFT ] ) )
                        lstr_tab_dependencies.append( self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
                        lstr_tab_dependencies.append( self.func_convert_fastq_rna_bam( dict_sample[ STR_RNA_LEFT ] ) )
                        lstr_tab_dependencies.append( self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                        lstr_generated_dna_vcf_snps.append( self.func_convert_fastq_dna_vcf( dict_sample[ STR_DNA_LEFT ] ) )
                        lstr_generated_rna_vcf_snps.append( self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) )
                        lstr_generated_dna_rna_tab.append( self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                        llstr_tab_samples.append( [ dict_sample[ STR_MAF_SAMPLE_NAME ],
                                                    "NA",
                                                    "NA",
                                                    self.func_convert_fastq_rna_bam( dict_sample[ STR_RNA_LEFT ] ),
                                                    self.func_convert_fastq_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) + "\n" ] )
 
                with open( str_tab_sample_file, "w" ) as hndl_tab_samples:
                    hndl_tab_samples.write( llstr_tab_samples )

                # Make tabs for biological truth data and RNASEQ Data
                # Formats: DNASEQ only
                # ID from RNA Left \t DNA bam \t DNA vcf \t RNA bam \t RNA vcf
                # Formats: MAF only
                # MAF ID \t NA \t NA \t RNA bam \t RNA vcf
                # Formats: MAf and DNASEQ
                # MAF ID \t DNA bam \t DNA vcf \t RNA bam \t RNA vcf
                str_cmd_make_tabs = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annotation_config,
                                                "--run_conf", str_tab_run_conf_file, "--reads_list_file", str_tab_sample_file,
                                                "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parrallel_procs", arg_parsed.i_jobs ] )
                lcmd_commands.append( Command.Command( str_cur_command = str_cmd_make_tab,
                                                       lstr_cur_dependencies = lstr_tab_dependencies,
                                                       lstr_products = lstr_tab_products ) )

            # Make synthetic comparisons tab files
            if dict_sample[ STR_SYNTHETIC_BAM ]:
                # Make tabs for synthetic truth data and RNASEQ Data
                # Formats: Synthetic
                # ID from RNA left \t Synthetic bam \t Synthetic vcf \t RNA bam \t RNA vcf
                str_synthetic_tab_sample_file = os.path.join( str_current_project_dir, "tab_syn_samples.txt" )
                lstr_syn_dependencies = []
                lstr_syn_products = []
                for dict_sample in dict_sample_study.items():
                    if not dict_sample[ STR_SYNTHETIC_BAM ] or not dict_sample[ STR_SYNTHETIC_VCF ]:
                        continue
                    llstr_synthetic_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0],
                                              dict_sample[ STR_SYNTHETIC_BAM ],
                                              dict_sample[ STR_SYNTHETIC_VCF ],
                                              self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ] ),
                                              self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) + "\n" ] ]
                    lstr_syn_tab_dependencies.extend( [ dict_sample[ STR_SYNTHETIC_BAM ], dict_sample[ STR_SYNTHETIC_VCF ],
                                                 self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ] ),
                                                 self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ] ) ] )
                    lstr_syn_tab_products.append( func_convert_fastq_syn_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                    lstr_generated_syn_rna_tab.append(  func_convert_fastq_syn_rna_tab( dict_sample[ STR_RNA_LEFT ] ) )
                with open( str_synthetic_tab_sample_file, "w" ) as hndl_syn_samples:
                    hndl_syn_samples.write( llstr_synthetic_samples )
                str_cmd_make_synthetic_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annotation_config,
                                                "--run_conf", str_synthetic_tab_run_conf_file, "--reads_list_file", str_synthetic_tab_sample_file,
                                                "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parrallel_procs", arg_parsed.i_jobs ] )
                lcmd_commands.append( Command.Command( str_cur_command = str_cmd_make_synthetic_tab,
                                                       lstr_cur_dependencies = lstr_syn_tab_dependencies,
                                                       lstr_products = lstr_syn_tab_products ) )

            # Make figures for validation with a study
            dict_validation_commands = func_validation_figure_commands( args_parsed = args_parsed, str_cur_project_dir = str_cur_project_dir,
                                                                        lstr_dna_vcf_snps = lstr_generated_vcf_snps,
                                                                        lstr_rna_vcfs_snp = lstr_generated_vcf_snps,
                                                                        lstr_maf_rna_tab = lstr_generated_maf_rna_tab,
                                                                        lstr_maf_dna_tab = lstr_generated_maf_dna_tab,
                                                                        lstr_dna_rna_tab = lstr_generated_dna_rna_tab,
                                                                        lstr_syn_rna_tab = lstr_generated_syn_rna_tab ) 
            lcmd_commands.extend( dict_validation_commands[ STR_CMDS ] )

        # TODO Add in figures for validation between studies    
        # This includes ROCs of methods side by side so ROCs for N! number of run confs.
        return( lcmd_commands )

    def func_validation_figure_commands( self, args_parsed, str_cur_project_dir, lstr_dna_vcf_snps = [], lstr_rna_vcfs_snp = [], lstr_maf_rna_tab = [], lstr_maf_dna_tab = [], lstr_dna_rna_tab = [], lstr_syn_rna_tab = [] ):
        """

        * return : List of commands to be ran
                 : list of command objects
        """

        # Collects commands
        lcmd_commands = []

        # Distance metric functions
        str_jaccard_distance = os.path.join( self.str_src_dir, "jaccard_distance.R" )
        str_jaccard_distance_restricted = os.path.join( self.str_src_dir, "jaccard_distance_restricted.R" )

        # Check to make sure directories exit
        STR_FIGURE_DIR = os.path.join( str_cur_project_dir, self.str_figure_dir )
        STR_DNA_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_dna )
        STR_RNA_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_rna )
        STR_MAF_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dir, self.str_figure_maf )

        # Plot the rate of mutation over a window.
        # Mutation rates are often local.
        # Would be interesting to filter on # of mutations in a gene given the mutation rates locally.
        if( len( lstr_dna_vcfs_snps ) > 1 and len( lstr_rna_vcfs_snps ) > 1 ):
            # Create genotype matrix and list
            ## DNA and RNA
            str_genotype_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype_matrix.txt" )
            str_cmd_create_genotype_matrix = os.path.join( self.str_src_dir, "vcfs_to_genotype_matrix.py" ) + " --matrix " + str_genotype_matrix + " " + " ".join( lstr_rna_vcfs_snps + lstr_dna_vcfs_snps )
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_create_genotype_matrix,
                                               lstr_cur_dependencies = lstr_rna_vcfs_snps + lstr_dna_vcfs_snps,
                                               lstr_cur_products = [ str_genotype_matrix ] ) )
            # Visualize genotype matrices
            ## DNA and RNA, restricted by type
            str_genotype_matrix_restricted_pdf = os.path.join( STR_FIGURE_DIR, "sample_dna_rna_genotype_matrix_restricted.pdf" )
            str_genotype_distance_restricted_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype_restricted.dist" )
            str_cmd_genotype_matrix_dna_rna_restricted = " ".join( [ os.path.join( self.str_src_dir, "make_dendrogram_generic.R" ), "--input_matrix", str_genotype_matrix, 
                                                      "--distance_function", str_jaccard_distance_restricted,
                                                      "--output_pdf", str_genotype_matrix_restricted_pdf,
                                                      "--output_distance_matrix", str_genotype_distance_restricted_matrix ] ) 
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_dna_rna_restricted,
                                               lstr_cur_dependencies = [ str_genotype_matrix, str_jaccard_distance_restricted ],
                                               lstr_cur_products = [ str_genotype_matrix_restricted_pdf, str_genotype_distance_restricted_matrix ] ) )
            ## DNA and RNA not restricted by type
            str_genotype_matrix_pdf = os.path.join( STR_FIGURE_DIR, "sample_dna_rna_genotype_matrix.pdf" )
            str_genotype_distance_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype.dist" )
            str_cmd_genotype_matrix_dna_rna = " ".join( [ os.path.join( self.str_src_dir, "make_dendrogram_generic.R" ), "--input_matrix", str_genotype_matrix, 
                                                      "--distance_function", str_jaccard_distance,
                                                      "--output_pdf", str_genotype_matrix_pdf,
                                                      "--output_distance_matrix", str_genotype_distance_matrix ] ) 
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_dna_rna,
                                               lstr_cur_dependencies = [ str_genotype_matrix, str_jaccard_distance ],
                                               lstr_cur_products = [ str_genotype_matrix_pdf, str_genotype_distance_matrix ] ) )

        if( len( lstr_rna_vcfs_snps ) > 1 ):
            ## RNA
            str_genotype_rna_matrix = os.path.join( args_parsed.str_file_base, "sample_rna_genotype_matrix.txt" )
            str_cmd_create_genotype_rna_matrix = os.path.join( self.str_src_dir, "vcfs_to_genotype_matrix.py" ) + " --matrix " + str_genotype_rna_matrix + " " + " ".join( lstr_rna_vcfs_snps )
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_create_genotype_rna_matrix,
                                               lstr_cur_dependencies = lstr_rna_vcfs_snps,
                                               lstr_cur_products = [ str_genotype_rna_matrix ] ) )
        
            # Visualize genotype matrices
            ## RNA
            str_genotype_rna_matrix_pdf = os.path.join( STR_FIGURE_DIR, "sample_rna_genotype_matrix.pdf" )
            str_genotype_rna_distance_matrix = os.path.join( args_parsed.str_file_base, "sample_rna_genotype.dist" )
            str_cmd_genotype_matrix_rna = " ".join( [ os.path.join( self.str_src_dir, "make_dendrogram_generic.R" ), "--input_matrix", str_genotype_rna_matrix, 
                                                      "--distance_function", str_jaccard_distance,
                                                      "--output_pdf", str_genotype_rna_matrix_pdf,
                                                      "--output_distance_matrix", str_genotype_rna_distance_matrix ] ) 
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_rna,
                                               lstr_cur_dependencies = [ str_genotype_rna_matrix, str_jaccard_distance ],
                                               lstr_cur_products = [ str_genotype_rna_matrix_pdf, str_genotype_rna_distance_matrix ] ) )

        # Explore false positive and negative rates
        # DNA_RNA
        str_DNA_RNA_figures = os.path.join( STR_FIGURE_DIR, "dna_rna" )
        lstr_DNA_RNA_roc_1 = [ os.path.join( str_DNA_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_dna_rna_tab ]
        lstr_DNA_RNA_roc_2 = [ os.path.join( str_DNA_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_dna_rna_tab ]
        str_compare_DNA_RNA_cmd = os.path.join( self.str_src_dir, "visualize_mutation_depth_tab_files.R" ) + " -o " + str_DNA_RNA_figures + " -k DNA_RNA --serial_plots " + " ".join( lstr_dna_rna_tab )
        lcmd_commands.append( Command.Command( str_cur_command = str_compare_DNA_RNA_cmd,
                                               lstr_cur_dependencies = lstr_dna_rna_tab,
                                               lstr_cur_products = lstr_DNA_RNA_roc_1 + lstr_DNA_RNA_roc_2 ) )

        if not args_parsed.str_maf_file is None:
            # MAF_DNA
            str_MAF_DNA_figures = os.path.join( STR_FIGURE_DIR, "maf_dna" )
            lstr_MAF_DNA_roc_1 = [ os.path.join( str_MAF_DNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_maf_dna_tab ]
            lstr_MAF_DNA_roc_2 = [ os.path.join( str_MAF_DNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_maf_dna_tab ]
            str_compare_MAF_DNA_cmd = os.path.join( self.str_src_dir, "visualize_mutation_depth_tab_files.R" ) + " -o " + str_MAF_DNA_figures + " -k MAF_DNA --serial_plots " + " ".join( lstr_maf_dna_tab )
            lcmd_commands.append( Command.Command( str_cur_command = str_compare_MAF_DNA_cmd,
                                                   lstr_cur_dependencies = lstr_maf_dna_tab,
                                                   lstr_cur_products = lstr_MAF_DNA_roc_1 + lstr_MAF_DNA_roc_2 ) )

            # MAF_RNA
            str_MAF_RNA_figures = os.path.join( STR_FIGURE_DIR, "maf_rna" )
            lstr_MAF_RNA_roc_1 = [ os.path.join( str_MAF_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_10_pred_vary.pdf" ] ) ) for str_file in lstr_maf_rna_tab ]
            lstr_MAF_RNA_roc_2 = [ os.path.join( str_MAF_RNA_figures, "_".join( [ os.path.basename( str_file ), "roc_truth_vary_pred_1.pdf" ] ) ) for str_file in lstr_maf_rna_tab ]
            str_compare_MAF_RNA_cmd = os.path.join( self.str_src_dir, "visualize_mutation_depth_tab_files.R" ) + " -o " + str_MAF_RNA_figures + " -k MAF_RNA --serial_plots " + " ".join( lstr_maf_rna_tab )
            lcmd_commands.append( Command.Command( str_cur_command = str_compare_MAF_RNA_cmd,
                                                   lstr_cur_dependencies = lstr_maf_rna_tab,
                                                   lstr_cur_products = lstr_MAF_RNA_roc_1 + lstr_MAF_RNA_roc_2 ) )

            # MAF
            str_key_mutations_tab = "".join( [ " -t " + str_tab for str_tab in lstr_maf_dna_tab ] )
            str_key_mutation_MAF_DNA_output_pdf = os.path.join( STR_MAF_FIGURE, "key_mutation_counts.pdf" )
            str_key_mutations_cmd = os.path.join( self.str_src_dir, "tabs_to_percent_mutations.py" ) + " --gtf " + args_parsed.str_annotation_gtf + " --key " + args_parsed.str_key_mutations + " -o " + str_key_mutation_MAF_DNA_output_pdf + str_key_mutations_tab
            lcmd_commands.append( Command.Command( str_cur_command = str_key_mutations_cmd,
                                               lstr_cur_dependencies = lstr_maf_dna_tab + [ args_parsed.str_annotation_gtf ],
                                               lstr_cur_products = [ str_key_mutation_MAF_DNA_output_pdf ] ) )
            # RNA
            str_key_mutations_tab = "".join( [ " -t " + str_tab for str_tab in lstr_maf_rna_tab ] )
            str_key_mutation_MAF_RNA_output_pdf = os.path.join( STR_RNA_FIGURE, "key_mutation_counts.pdf" )
            str_key_mutations_cmd = os.path.join( self.str_src_dir, "tabs_to_percent_mutations.py" ) + " --second --gtf " + args_parsed.str_annotation_gtf + " --key " + args_parsed.str_key_mutations + " -o " + str_key_mutation_MAF_RNA_output_pdf + str_key_mutations_tab
            lcmd_commands.append( Command.Command( str_cur_command = str_key_mutations_cmd,
                                                   lstr_cur_dependencies = lstr_maf_rna_tab + [ args_parsed.str_annotation_gtf ],
                                                   lstr_cur_products = [ str_key_mutation_MAF_RNA_output_pdf ] ) )

        # Count mutations and plot
        # DNA
        str_key_mutations_tab = "".join( [ " -t " + str_tab for str_tab in lstr_dna_rna_tab ] )
        str_key_mutation_DNA_RNA_output_pdf = os.path.join( STR_DNA_FIGURE, "key_mutation_counts.pdf" )
        str_key_mutations_cmd = os.path.join( self.str_src_dir, "tabs_to_percent_mutations.py" ) + " --gtf " + args_parsed.str_annotation_gtf + " --key " + args_parsed.str_key_mutations + " -o " + str_key_mutation_DNA_RNA_output_pdf + str_key_mutations_tab
        lcmd_commands.append( Command.Command( str_cur_command = str_key_mutations_cmd,
                                               lstr_cur_dependencies = lstr_dna_rna_tab + [ args_parsed.str_annotation_gtf ],
                                               lstr_cur_products = [ str_key_mutation_DNA_RNA_output_pdf ] ) )

        return lcmd_commands

if __name__ == "__main__":

    # Run pipeline
    RNASEQ_mutation_validation().func_run_pipeline()
