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
import sciedpiper.PipelineRunner as PipelineRunner

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

# Keys for the JSON object that interfaces with the mutation isnpector visualization app
STR_INSPECTOR_TP = "TP"
STR_INSPECTOR_FP = "FP"
STR_INSPECTOR_FN = "FN"
STR_INSPECTOR_RNA = "RNA"
STR_INSPECTOR_DNA = "DNA"

# Paired sample file indices
I_SYNTHETIC_VCF = 5
I_RNA_RIGHT = 4
I_RNA_LEFT = 3
I_DNA_RIGHT = 2
I_DNA_LEFT = 1
I_MAF_SAMPLE = 0

class RNASEQ_mutation_validation(PipelineRunner.PipelineRunner):

    def __init__( self ):

        # File structure
        self.str_bams_dir = "alignment"
        self.str_conf = "conf"
        self.str_figure_dir = "figure"
        self.str_log_dir = "log"
        self.str_sample_files_dir = "sample_files"
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
        arg_raw.add_argument( "--align_conf", dest = "str_align_run_config", help = "Alignment config for all variant calling." )

        # Annotation config file for all runs
        arg_raw.add_argument( "--annot_conf", dest = "str_annot_config", required = True, help = "Annotation config file for all samples." )

        # Conf file, one run or all samples per conf file will be ran
        arg_raw.add_argument( "--call_conf", dest = "lstr_run_configs", action = "append", default = [], help = "Config file for the call variants part of the pipeline. Can be used multiple times, each conf file is a setting for a study that is ran on ALL samples and validated as a seperate group." )

        # Number of threads
        arg_raw.add_argument( "--jobs", dest = "i_jobs", default = 1, action = "store", help = "Max number of jobs to run at a time." )

        # Truth calling run conf happends once (on all DNA samples) per validation run.
        arg_raw.add_argument( "--truth_call_conf", dest = "str_truth_run_config", default = None, action = "store", help = "Config file to call variants on the truth data." )

        # File that pairs RNA and DNA seq files
        arg_raw.add_argument( "--sample_file", dest = "str_sample_pairing_file", required = True, help = "A file that contains pairing of RNA and DNA Bam and Vcf files." )

        # MAF file
        # dest = str_maf_file
        arg_raw.add_argument( "--maf", dest = "str_maf_file", default = None, help = "The maf file used to compare DNA dn RNA vcfs." )

        # Reduced dbSNP vcf file
        # dest = str_reduced_dbsnp_vcf
        arg_raw.add_argument( "--ref_vcf", dest = "str_reduced_dbsnp_vcf", help = "A reference vcf file." )
        arg_raw.add_argument("--ref_gtf", dest = "str_annotation_gtf", help = "A reference gtf file." )

        # Comma delimited string of key mutations to look at occurence in samples
        # dest = args_parsed.str_key_mutations
        arg_raw.add_argument( "--key", dest = "str_key_mutations", help = "Comma delimited string of key mutations (gene names) to look at occurence in samples." )

        # Run config file for making tab files
        arg_raw.add_argument( "--tab_conf_dir", dest = "str_tab_config_dir", action = "store", help = "Directory to find run config file for making tab files" )

    def func_convert_fastq_left_rna_bam( self, str_fastq_left, str_output_dir ):
        # Convert left fastq rna sample to the commonly aligned bam (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
#        return os.path.join( str_output_dir, self.str_bams_dir, "samples", str_file_base, "COMMON_ALIGNMENT", "star_align_2_reads_"+os.path.basename( str_fastq_left ).replace(".","_")+"_P_qtrim_fq","Aligned.sorted.bam" )
        return os.path.join( str_output_dir, self.str_bams_dir, "samples", str_file_base.replace(".","_"), "COMMON_ALIGNMENT", "misc","Aligned.sorted.bam" )

    def func_convert_fastq_left_rna_vcf( self, str_fastq_left, str_output_dir ):
        # Convert left fastq rna sample to the final annotated vcf (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "RNASEQ_MUTATION_PREMADE", "variants_annotated.vcf.gz" )

    def func_convert_fastq_left_cancer_tab( self, str_fastq_left, str_output_dir ):
        # Convert left fastq rna sample to the final cancer tab (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "RNASEQ_MUTATION_PREMADE", "cancer.tab" )

    def func_convert_fastq_left_cancer_vcf( self, str_fastq_left, str_output_dir ):
        # Convert left fastq rna sample to the final cancer vcf (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "RNASEQ_MUTATION_PREMADE", "cancer.vcf" )

    def func_convert_fastq_left_dna_bam( self, str_fastq_left, str_output_dir ):
        # Convert left fastq dna sample name to the commonly aligned bam
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
#        return os.path.join( str_output_dir, "truth_runs", "samples", str_file_base, "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", "misc_reads_"+str_file_base_misc+"_fq_P_qtrim_fq", "reads_" + str_file_base_misc + ".recal_snp.bam" )
        return os.path.join( str_output_dir, "truth_runs", "samples", str_file_base.replace(".","_"), "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", "misc", "reads_" + str_file_base.replace(".","_") + ".recal_snp.bam" )

    def func_convert_fastq_left_dna_vcf( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq dna sample to the fnal output vcf
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
#        return os.path.join( str_output_dir, "truth_runs", "samples", str_file_base, "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", "misc_reads_"+str_file_base_misc+"_fq_P_qtrim_fq", "variants_annotated.vcf" )
        return os.path.join( str_output_dir, "truth_runs", "samples", str_file_base.replace(".","_"), "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", "misc", "reads_"+str_file_base.replace(".","_")+".variants_filtered.vcf" )

    def func_convert_fastq_maf_dna_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq dna to the maf vs dna tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_MAF_DNA", "maf_dna.tab" )

    def func_convert_fastq_maf_rna_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq rna to the maf vs rna tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_MAF_RNA", "maf_rna.tab" )

    def func_convert_fastq_dna_rna_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq dna to the dna vs rna tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_DNA_RNA", "dna_rna.tab" )

    def func_convert_fastq_dna_rna_edit_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq dna to the dna vs rna tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_DNA_RNA_EDIT", "dna_rna_rnaedit.tab" )

    def func_convert_fastq_dna_rna_cosmic_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq dna to the dna vs rna tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_DNA_RNA_COSMIC", "dna_rna_cosmic.tab" )

    def func_convert_fastq_dna_rna_cancer_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq dna to the dna vs rna tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_DNA_RNA_CANCER", "dna_rna_cancer.tab" )

    def func_convert_fastq_syn_rna_tab( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq rna sample to the synthetic tab sample (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "TABS_SYN", "syn.tab" )

    def func_convert_fastq_rna_depth_gatk( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq rna sample to the depth file for the common alignment bam (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, "samples", str_file_base.replace(".","_"), "RNASEQ_MUTATION_PREMADE", "misc", str_file_base + ".depth" )

    def func_convert_fastq_rna_depth_samtools( self, str_fastq_left, str_output_dir ):
        # Convert the left fastq rna sample to the depth file for the common alignment bam (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, self.str_bams_dir, "samples", str_file_base.replace(".","_"), "COMMON_ALIGNMENT", "Aligned.depth" )

    def func_convert_fastq_dna_depth_gatk( self, str_fastq_left, str_output_dir ):
        # Convert the right fastq dna sample to the depth file fo the common alignment bam (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, self.str_truth_runs, "samples", str_file_base.replace(".","_"), "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", str_file_base + ".depth" )

    def func_convert_fastq_dna_depth_samtools( self, str_fastq_left, str_output_dir ):
        # Convert the right fastq dna sample to the depth file fo the common alignment bam (ok)
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_output_dir, self.str_truth_runs, "samples", str_file_base.replace(".","_"), "COMMON_ALIGNMENT", "Aligned.depth" )

    def func_convert_fastq_left_study_dir( self, str_fastq_left, str_parent_dir ):
        str_file_base = os.path.splitext( os.path.basename( str_fastq_left ) )[0]
        return os.path.join( str_parent_dir, str_file_base.replace(".","_") )

    def func_run_truth_calling( self, args_parsed, dict_sample_study ):

        lstr_truth_call_products = []
        lcmd_truth_calling = []

        # Make DNA Seq / truth samples.txt
        str_truth_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "truth_samples.txt" )
        llstr_truth_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_DNA_LEFT ] ) )[0].replace(".","_"),
                                  dict_sample[ STR_DNA_LEFT],
                                  dict_sample[ STR_DNA_RIGHT ],
                                  self.func_convert_fastq_left_study_dir( dict_sample[ STR_DNA_LEFT ], os.path.join( args_parsed.str_file_base, self.str_truth_runs, "samples" ) ) ]
            for dict_sample in dict_sample_study.values() if dict_sample[ STR_DNA_LEFT ] and dict_sample[ STR_DNA_RIGHT ] ]
        lstr_truth_samples = [ "\t".join( lstr_files ) for lstr_files in llstr_truth_samples ]
        if len( lstr_truth_samples ) > 0:
            if not args_parsed.f_Test:
                with open( str_truth_sample_file, "w" ) as hndl_truth_samples:
                    hndl_truth_samples.write( "\n".join( lstr_truth_samples ) )

            # Make dep / product names
            lstr_truth_call_dependencies = [ [ lstr_sample[ 1 ], lstr_sample[2] ] for lstr_sample in llstr_truth_samples ]
            lstr_truth_call_dependencies = [ str_dep for lstr_subgroup in lstr_truth_call_dependencies for str_dep in lstr_subgroup ]
            lstr_truth_call_products = [ [ self.func_convert_fastq_left_dna_bam( lstr_file[ 1 ], args_parsed.str_file_base ),
                                           os.path.splitext( self.func_convert_fastq_left_dna_bam( lstr_file[ 1 ], args_parsed.str_file_base ) )[0] + ".bai",
                                           self.func_convert_fastq_left_dna_vcf ( lstr_file[ 1 ], args_parsed.str_file_base )] for lstr_file in llstr_truth_samples ]
            lstr_truth_call_products = [ str_prod for lstr_subgroup in lstr_truth_call_products for str_prod in lstr_subgroup ]

            # Run dna samples
            str_cmd_truth_call_variants = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                            "--run_conf", args_parsed.str_truth_run_config, "--reads_list_file", str_truth_sample_file, 
                                            "--project_base_dir", os.path.join( args_parsed.str_file_base, self.str_truth_runs ),
                                            "--memory 35 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                            #"--memory 35 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
            # TODO lcmd_truth_calling.append( Command.Command( str_cur_command = str_cmd_truth_call_variants,
            # TODO                                   lstr_cur_dependencies = lstr_truth_call_dependencies + [ args_parsed.str_annot_config, args_parsed.str_truth_run_config ],
            # TODO                                   lstr_cur_products = lstr_truth_call_products ) )
        return lcmd_truth_calling

    def func_prep_reference_vcf( self, args_parsed ):
        str_filtered_dbsnp_vcf = os.path.join( args_parsed.str_file_base, self.str_filtered_vcf, os.path.splitext( os.path.basename( args_parsed.str_reduced_dbsnp_vcf ) )[ 0 ] + "_snp.vcf" )
        str_filtered_dbsnp_vcf_command = os.path.join( "reduce_vcf_to_snps.py --reference " ) + args_parsed.str_reduced_dbsnp_vcf + " " + str_filtered_dbsnp_vcf
        return [ Command.Command( str_cur_command = str_filtered_dbsnp_vcf_command,
                                               lstr_cur_dependencies = [ args_parsed.str_reduced_dbsnp_vcf ],
                                               lstr_cur_products = [ str_filtered_dbsnp_vcf ] ) ]

    def func_make_snp_calling_sample_file( self, args_parsed, dict_sample_study, str_project_dir ):
        str_variant_calling_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, os.path.basename( str_project_dir ) + ".txt" )
        llstr_bam_call_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0].replace(".","_"),
                                     self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT], args_parsed.str_file_base),
                                     self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT], args_parsed.str_file_base ),
                                     os.path.join( str_project_dir, "samples", os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0].replace(".","_") ) ]
                                  for dict_sample in dict_sample_study.values() if dict_sample[ STR_RNA_LEFT ] ]
        lstr_bam_call_samples = [ "\t".join( lstr_files ) for lstr_files in llstr_bam_call_samples ]
        # Make variant calling sample file
        if len( lstr_bam_call_samples ) > 0 and not args_parsed.f_Test:
            with open( str_variant_calling_sample_file, "w" ) as hndl_bam_samples:
                hndl_bam_samples.write( "\n".join( lstr_bam_call_samples ) )
        return [ llstr_bam_call_samples, str_variant_calling_sample_file ]

    def func_make_common_alignment_sample_file( self, args_parsed, dict_sample_study ):
        str_align_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, "variant_calling_samples.txt" )
        llstr_align_samples = [ [ os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0].replace(".","_"),
                                  dict_sample[ STR_RNA_LEFT],
                                  dict_sample[ STR_RNA_RIGHT ],
                                  self.func_convert_fastq_left_study_dir( dict_sample[ STR_RNA_LEFT ], os.path.join( args_parsed.str_file_base, self.str_bams_dir, "samples" ) ) ]
                              for dict_sample in dict_sample_study.values() if dict_sample[ STR_RNA_LEFT ] and dict_sample[ STR_RNA_RIGHT ] ]
        lstr_align_samples = [ "\t".join( lstr_files ) for lstr_files in llstr_align_samples ]

        # Make rnaseq alignment sample file
        if len( lstr_align_samples ) > 0 and not args_parsed.f_Test:
            with open( str_align_sample_file, "w" ) as hndl_rna_samples:
                hndl_rna_samples.write( "\n".join( lstr_align_samples ) )
        return [ llstr_align_samples, str_align_sample_file ]

    def func_make_link_command( self, str_original_file, str_new_file ):
        str_command = "ln -sf " + str_original_file + " " + str_new_file
        cmd_ln =  Command.Command( str_cur_command = str_command,
                                                   lstr_cur_dependencies = [ str_original_file ],
                                                   lstr_cur_products = [ str_new_file ] )
        return cmd_ln

    ############## These will change, not working with a proper dataset
    def func_make_cancer_filtered_vcf_path( self, str_sample_dir, str_sample_name ):
        return os.path.join( str_sample_dir, "RNASEQ_MUTATION_PREMADE", "cancer.vcf"  )

    def func_make_cosmic_filtered_vcf_path( self, str_sample_dir, str_sample_prefix ):
        return os.path.join( str_sample_dir, "RNASEQ_MUTATION_PREMADE", "misc", str_sample_prefix + "_clean_snp_RNAedit.vcf_dbsnp_cosmic_filtered.vcf" )

    def func_make_rnaediting_filtered_vcf_path( self, str_sample_dir, str_sample_prefix ):
        return os.path.join( str_sample_dir, "RNASEQ_MUTATION_PREMADE", "misc", str_sample_prefix + "_clean_snp_RNAedit.vcf" )

    def func_make_initial_filtered_vcf_path( self, str_sample_dir, str_sample_prefix ):
        return os.path.join( str_sample_dir, "RNASEQ_MUTATION_PREMADE", "misc", str_sample_prefix + "_clean_snp.vcf" )

    def func_make_calling_depth_path( self, str_sample_dir ):
        return os.path.join( str_sample_dir, "RNASEQ_MUTATION_PREMADE", "RNASEQ_MUTATION_PREMADE.depth" )

    def func_make_truth_vcf_path( self, str_sample_dir, str_sample_name ):
        return os.path.join( str_sample_dir, "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", "misc", "reads_"+str_sample_name+".variants_filtered.vcf" )

    def func_make_truth_depth_path( self, str_sample_dir ):
        return os.path.join( str_sample_dir, "PREMADE_GATK_DNASEQ_LIKE_RNASEQ", "PREMADE_GATK_DNASEQ_LIKE_RNASEQ.depth" )
    #################################################################################

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
                                          for str_file in [ self.str_bams_dir, self.str_log_dir, 
                                                            self.str_sample_files_dir, self.str_truth_runs, self.str_filtered_vcf ]] ): 
            exit( 2 )

        #########################################
        # Prep secondary files if needed.
        # SOFT Filter the DBSNP vcf file to SNP
        #########################################
        lcmd_commands_run.extend( self.func_prep_reference_vcf( args_parsed ) )

        ###############################
        # Run the truth runs if needed
        ###############################
        # TODO if args_parsed.str_truth_run_config:
        # TODO     lcmd_commands_run.extend( self.func_run_truth_calling( args_parsed, dict_sample_study ) )

        ############################################################
        # Make sample file for common alignment and snp calling
        ############################################################
        llstr_alignment_samples, str_alignment_sample_file = self.func_make_common_alignment_sample_file( args_parsed, dict_sample_study )

        if args_parsed.lstr_run_configs and args_parsed.str_align_run_config:
            # Make dep / product names
            lstr_bam_dependencies = [ [ lstr_sample[ 1 ], lstr_sample[2] ] for lstr_sample in llstr_alignment_samples ]
            lstr_bam_dependencies = [ str_dep for lstr_subgroup in lstr_bam_dependencies for str_dep in lstr_subgroup ]
            lstr_bam_products = [ self.func_convert_fastq_left_rna_bam( lstr_sample[ 1 ], args_parsed.str_file_base ) for lstr_sample in llstr_alignment_samples ]
            lstr_bam_products.extend( [ str_bam + ".bai" for str_bam in lstr_bam_products ] )
            lstr_call_dependencies = lstr_bam_products

            ########################
            # Align RNA Seq bams
            ########################
            str_cmd_make_bam = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                            "--run_conf", args_parsed.str_align_run_config, "--reads_list_file", str_alignment_sample_file,
                                            "--project_base_dir", os.path.join( args_parsed.str_file_base, self.str_bams_dir),
                                            #"--memory 35 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                            "--memory 35 --num_parallel_procs", str( args_parsed.i_jobs ) ] )

            # TODO lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_make_bam,
            # TODO                                        lstr_cur_dependencies = lstr_bam_dependencies + [ args_parsed.str_annot_config, args_parsed.str_align_run_config ],
            # TODO                                        lstr_cur_products = lstr_bam_products ) )

            # For each conf file
            for str_call_run_conf in args_parsed.lstr_run_configs:
                # Current project directory
                str_call_basename = os.path.splitext( os.path.basename( str_call_run_conf ) )[ 0 ]
                str_current_project_dir = os.path.join( args_parsed.str_file_base, str_call_basename )
                if not cur_pipeline.func_mkdirs( [ str_current_project_dir ] ):
                    exit( 3 )
                # Create sample file for the calling
                llstr_calling_samples, str_calling_sample_file = self.func_make_snp_calling_sample_file( args_parsed, dict_sample_study, str_current_project_dir )
                if len( [ lstr_calling_samples for lstr_calling_samples in llstr_calling_samples if lstr_calling_samples ] ) < 1:
                    continue
                # Track files made
                lstr_generated_dna_vcf_snps = []
                lstr_generated_rna_vcf_snps = []
                lstr_generated_maf_rna_tab = []
                lstr_generated_maf_dna_tab = []
                lstr_generated_dna_rna_tab = []
                lstr_generated_dna_rna_tab_cancer = []
                lstr_generated_dna_rna_tab_cosmic = []
                lstr_generated_dna_rna_tab_edit = []
                lstr_generated_syn_rna_tab = []

                #############################
                # Call variants off of a bam
                #############################
                str_cmd_call_variants = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", str_call_run_conf, "--reads_list_file", str_calling_sample_file, 
                                                "--project_base_dir", str_current_project_dir, "--memory 35 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                                #"--project_base_dir", str_current_project_dir, "--memory 35 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                lstr_call_products = [ self.func_convert_fastq_left_rna_vcf( lstr_sample[ 1 ], str_current_project_dir ) for lstr_sample in llstr_alignment_samples ]
                lstr_call_products = lstr_call_products + [ self.func_convert_fastq_left_cancer_vcf( lstr_sample[ 1 ], str_current_project_dir ) for lstr_sample in llstr_alignment_samples ]
                lstr_call_products = lstr_call_products + [ self.func_convert_fastq_left_cancer_tab( lstr_sample[ 1 ], str_current_project_dir ) for lstr_sample in llstr_alignment_samples ]

                # TODO lcmd_commands_run.append( Command.Command( str_cur_command = str_cmd_call_variants,
                # TODO                                    lstr_cur_dependencies = lstr_call_dependencies + [ args_parsed.str_annot_config, str_call_run_conf ],
                # TODO                                    lstr_cur_products = lstr_call_products ) )

                ##########################################
                # Make dna rna comparisons tab files
                # Make tabs for Exome seq and RNASEQ Data
                ##########################################
                if args_parsed.str_tab_config_dir:
                    # Sample file
                    str_dna_rna_tab_sample_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, str_call_basename + "_tab_dna_rna_samples.txt" )
                    llstr_dna_rna_samples = []
                    # Config files for each ROC comparison
                    str_dna_rna_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, os.path.basename( os.path.splitext( str_call_run_conf )[0] ) + "_dna_rna_tabs.conf" )
                    str_dna_rna_edit_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, os.path.basename( os.path.splitext( str_call_run_conf )[0] ) + "_dna_rna_edit_tabs.conf" )
                    str_dna_rna_cosmic_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, os.path.basename( os.path.splitext( str_call_run_conf )[0] ) + "_dna_rna_cosmic_tabs.conf" )
                    str_dna_rna_cravat_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, os.path.basename( os.path.splitext( str_call_run_conf )[0] ) + "_dna_rna_cravat_tabs.conf" )
                    # Prod / Dep for str_dna_rna_tab_run_conf_file
                    lstr_dna_rna_tab_products = []
                    lstr_dna_rna_tab_edit_products = []
                    lstr_dna_rna_tab_cosmic_products = []
                    lstr_dna_rna_tab_cancer_products = []
                    lstr_generic_tab_dependencies = []
                    lstr_dna_rna_dependencies = []
                    lstr_dna_rna_edit_dependencies = []
                    lstr_dna_rna_cosmic_dependencies = []
                    lstr_dna_rna_cravat_dependencies = []

                    # Maf mapping info to make maf mapping files between vcf and maf keyword names
                    lstr_maf_mapping_init = []
                    lstr_maf_mapping_edit = []
                    lstr_maf_mapping_cosmic = []
                    lstr_maf_mapping_cancer = []
                    str_maf_mapping_init_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, str_call_basename + "_maf_init.mapping_txt" )
                    str_maf_mapping_edit_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, str_call_basename + "_maf_edit.mapping_txt" )
                    str_maf_mapping_cosmic_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, str_call_basename + "_maf_cosmic.mapping_txt" )
                    str_maf_mapping_cancer_file = os.path.join( args_parsed.str_file_base, self.str_sample_files_dir, str_call_basename + "_maf_cancer.mappig_txt" )

                    # JSON file prod/dep
                    lstr_dna_rna_file_info = []
                    lstr_dna_rna_json_file_dependencies = []

                    for dict_sample in dict_sample_study.values():
                        if not dict_sample[ STR_DNA_LEFT ]:
                            continue
                        str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0].replace(".","_")
                        str_cur_ref_sample = os.path.splitext( os.path.basename( dict_sample[ STR_DNA_LEFT ] ) )[0].replace(".","_")
                        str_calling_sample_dir = os.path.join( str_current_project_dir, "samples", str_cur_tab_sample )
                        str_truth_sample_dir = self.func_convert_fastq_left_study_dir( dict_sample[ STR_DNA_LEFT ], os.path.join( args_parsed.str_file_base, self.str_truth_runs, "samples" ) )
                        # Base for the links that are made
                        str_link_base = os.path.join( str_current_project_dir, "samples", str_cur_tab_sample, "reads." + str_cur_tab_sample )
                        # Tab file being made
                        str_cur_tab_file = self.func_convert_fastq_dna_rna_tab( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                        str_cur_edit_tab_file = self.func_convert_fastq_dna_rna_edit_tab( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                        str_cur_cosmic_tab_file = self.func_convert_fastq_dna_rna_cosmic_tab( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                        str_cur_cancer_tab_file = self.func_convert_fastq_dna_rna_cancer_tab( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                        llstr_dna_rna_samples.append( "\t".join( [ str_cur_tab_sample,
                                                                 os.path.join( str_current_project_dir, str_cur_tab_sample ), 
                                                                 os.path.join( str_current_project_dir, str_cur_tab_sample ) + "\n" ] ) )

                        # Make the sample prefix which differs between gatk and samtools
                        str_sample_link_prefix = "variants_filtered"
                        str_calling_link_prefix = "gatk"
                        if "samtools" in str_dna_rna_tab_run_conf_file.lower():
                            str_sample_link_prefix = "reads_"+str_cur_tab_sample+"_cluster_variants"
                            str_calling_link_prefix = "samtools"

                        # Links
                        # DNA reference depth
                        str_original_file = self.func_make_truth_depth_path( str_truth_sample_dir )
                        str_new_file = str_link_base + "_reference_dna_rna_"+str_calling_link_prefix+".depth"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_new_file ) )
                        lstr_generic_tab_dependencies.append( str_new_file )

                        # RNA depth
                        str_original_file = self.func_make_calling_depth_path( str_calling_sample_dir )
                        str_new_file = str_link_base + "_dna_rna_"+str_calling_link_prefix+".depth"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_new_file ) )
                        lstr_generic_tab_dependencies.append( str_new_file )

                        # DNA reference VCF
                        str_original_file = self.func_make_truth_vcf_path( str_truth_sample_dir, str_cur_ref_sample )
                        str_dna_rna_link_prod_ref_vcf = str_link_base + "_reference_dna_rna.vcf"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_dna_rna_link_prod_ref_vcf ) )
                        lstr_generic_tab_dependencies.append( str_dna_rna_link_prod_ref_vcf )
                        lstr_generated_dna_vcf_snps.append( str_dna_rna_link_prod_ref_vcf )

                        # RNA initial VCF
                        str_original_file = self.func_make_initial_filtered_vcf_path( str_calling_sample_dir, str_sample_link_prefix )
                        str_dna_rna_link_prod_vcf = str_link_base + "_" + str_calling_link_prefix + "_dna_rna.vcf"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_dna_rna_link_prod_vcf ) )
                        lstr_generated_rna_vcf_snps.append( str_dna_rna_link_prod_vcf )
                        lstr_dna_rna_dependencies.append( str_dna_rna_link_prod_vcf )
                        lstr_maf_mapping_init.append( [ dict_sample[ STR_MAF_SAMPLE_NAME ], str_original_file ] )

                        # RNA SNP Editing filter VCF
                        str_original_file = self.func_make_rnaediting_filtered_vcf_path( str_calling_sample_dir, str_sample_link_prefix )
                        str_new_file = str_link_base + "_" + str_calling_link_prefix + "_dna_rna_rnaedit.vcf"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_new_file ) )
                        lstr_dna_rna_edit_dependencies.append( str_new_file )
                        lstr_maf_mapping_edit.append( [ dict_sample[ STR_MAF_SAMPLE_NAME ], str_original_file ] )

                        # RNA COSMIC filter VCF
                        str_original_file = self.func_make_cosmic_filtered_vcf_path( str_calling_sample_dir, str_sample_link_prefix )
                        str_new_file = str_link_base + "_" + str_calling_link_prefix + "_dna_rna_cosmic.vcf"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_new_file ) )
                        lstr_dna_rna_cosmic_dependencies.append( str_new_file )
                        lstr_maf_mapping_cosmic.append( [ dict_sample[ STR_MAF_SAMPLE_NAME ], str_original_file ] )

                        # RNA Cancer filter VCF
                        str_original_file = self.func_make_cancer_filtered_vcf_path( str_calling_sample_dir, str_cur_tab_sample )
                        str_new_file = str_link_base + "_" + str_calling_link_prefix + "_dna_rna_cancer.vcf"
                        lcmd_commands_run.append( self.func_make_link_command( str_original_file, str_new_file ) )
                        lstr_dna_rna_cravat_dependencies.append( str_new_file )
                        lstr_maf_mapping_cancer.append( [ dict_sample[ STR_MAF_SAMPLE_NAME ], str_original_file ] )
                        
                        # Tab dependencies
                        lstr_dna_rna_tab_products.append( str_cur_tab_file )
                        lstr_dna_rna_tab_edit_products.append(str_cur_edit_tab_file )
                        lstr_dna_rna_tab_cosmic_products.append( str_cur_cosmic_tab_file )
                        lstr_dna_rna_tab_cancer_products.append(str_cur_cancer_tab_file )

                        # Tabs that will be made
                        lstr_generated_dna_rna_tab.append( str_cur_tab_file )
                        lstr_generated_dna_rna_tab_edit.append( str_cur_edit_tab_file )
                        lstr_generated_dna_rna_tab_cosmic.append( str_cur_cosmic_tab_file )
                        lstr_generated_dna_rna_tab_cancer.append( str_cur_cancer_tab_file )

                        # Add file info for the inspector's json object
                        # test_sample,test_rna.bam,test_dna.bam,test_rna.vcf,test_dna.vcf,test_comparison.tab
                        lstr_dna_rna_file_info.append( ",".join( [ str_cur_tab_sample, 
                                                           self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ),
                                                           self.func_convert_fastq_left_dna_bam( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ),
                                                           str_dna_rna_link_prod_vcf, str_dna_rna_link_prod_ref_vcf, str_cur_tab_file ] ) )
                        lstr_dna_rna_json_file_dependencies.extend( [ self.func_convert_fastq_left_rna_bam( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base ),
                                                           self.func_convert_fastq_left_dna_bam( dict_sample[ STR_DNA_LEFT ], args_parsed.str_file_base ),
                                                           str_dna_rna_link_prod_vcf, str_dna_rna_link_prod_ref_vcf, str_cur_tab_file ] )

                    # If there are samples to run make tab files.
                    if len( llstr_dna_rna_samples ) > 0:
                        if not args_parsed.f_Test:
                            with open( str_dna_rna_tab_sample_file, "w" ) as hndl_dna_rna_samples:
                                hndl_dna_rna_samples.write( "\n".join( llstr_dna_rna_samples ) )
                        ############### Make tab files
                        # Initial tab
                        str_cmd_make_dna_rna_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                    "--run_conf", str_dna_rna_tab_run_conf_file, "--reads_list_file", str_dna_rna_tab_sample_file,
                                                    #"--project_base_dir", str_current_project_dir, "--memory 20 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                                    "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                        #cmd_dna_rna_tab = Command.Command( str_cur_command = str_cmd_make_dna_rna_tab,
                        #                                   lstr_cur_dependencies = lstr_generic_tab_dependencies + lstr_dna_rna_dependencies + [ args_parsed.str_annot_config, str_dna_rna_tab_run_conf_file ],
                        #                                   lstr_cur_products = lstr_dna_rna_tab_products )
                        #lcmd_commands_run.append( cmd_dna_rna_tab )
                        # RNAEdit tab
                        str_cmd_make_dna_rna_edit_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                    "--run_conf", str_dna_rna_edit_tab_run_conf_file, "--reads_list_file", str_dna_rna_tab_sample_file,
                                                    #"--project_base_dir", str_current_project_dir, "--memory 20 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                                    "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                        #cmd_dna_rna_edit_tab = Command.Command( str_cur_command = str_cmd_make_dna_rna_edit_tab,
                        #                                   lstr_cur_dependencies = lstr_generic_tab_dependencies + lstr_dna_rna_edit_dependencies + [ args_parsed.str_annot_config, str_dna_rna_edit_tab_run_conf_file ],
                        #                                   lstr_cur_products = lstr_dna_rna_tab_edit_products )
                        #lcmd_commands_run.append( cmd_dna_rna_edit_tab )
                        # COSMIC tab
                        str_cmd_make_dna_rna_cosmic_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                    "--run_conf", str_dna_rna_cosmic_tab_run_conf_file, "--reads_list_file", str_dna_rna_tab_sample_file,
                                                    #"--project_base_dir", str_current_project_dir, "--memory 20 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                                    "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                        #cmd_dna_rna_cosmic_tab = Command.Command( str_cur_command = str_cmd_make_dna_rna_cosmic_tab,
                        #                                   lstr_cur_dependencies = lstr_generic_tab_dependencies + lstr_dna_rna_cosmic_dependencies + [ args_parsed.str_annot_config, str_dna_rna_cosmic_tab_run_conf_file ],
                        #                                   lstr_cur_products = lstr_dna_rna_tab_cosmic_products )
                        #lcmd_commands_run.append( cmd_dna_rna_cosmic_tab )
                        # Cancer tab
                        str_cmd_make_dna_rna_cravat_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                    "--run_conf", str_dna_rna_cravat_tab_run_conf_file, "--reads_list_file", str_dna_rna_tab_sample_file,
                                                    #"--project_base_dir", str_current_project_dir, "--memory 20 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                                    "--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                        #cmd_dna_rna_cravat_tab = Command.Command( str_cur_command = str_cmd_make_dna_rna_cravat_tab,
                        #                                   lstr_cur_dependencies = lstr_generic_tab_dependencies + lstr_dna_rna_cravat_dependencies + [ args_parsed.str_annot_config, str_dna_rna_cravat_tab_run_conf_file ],
                        #                                   lstr_cur_products = lstr_dna_rna_tab_cancer_products )
                        #lcmd_commands_run.append( cmd_dna_rna_cravat_tab )

                    ############################################
                    # Make figures for validation with a study
                    ############################################
                    if ( len( lstr_generated_dna_rna_tab ) + 
                         len( lstr_generated_dna_rna_tab_cosmic ) + 
                         len( lstr_generated_dna_rna_tab_cancer ) + 
                         len( lstr_generated_dna_rna_tab_edit ) +
                         len( lstr_dna_rna_dependencies ) +
                         len( lstr_dna_rna_edit_dependencies ) +
                         len( lstr_dna_rna_cosmic_dependencies ) +
                         len( lstr_dna_rna_cravat_dependencies ) ) > 0:

                        # If making figures, Write the maf mapping files to the project if workig with maf files.
                        if args_parsed.str_maf_file:
                            with open( str_maf_mapping_init_file, "w" ) as hndl_maf_init:
                                csv.writer( hndl_maf_init, delimiter="\t" ).writerows( lstr_maf_mapping_init )
                            with open( str_maf_mapping_edit_file, "w" ) as hndl_maf_edit:
                                csv.writer( hndl_maf_edit, delimiter="\t" ).writerows( lstr_maf_mapping_edit )
                            with open( str_maf_mapping_cosmic_file, "w" ) as hndl_maf_cosmic:
                                csv.writer( hndl_maf_cosmic, delimiter="\t" ).writerows( lstr_maf_mapping_cosmic )
                            with open( str_maf_mapping_cancer_file, "w" ) as hndl_maf_cancer:
                                csv.writer( hndl_maf_cancer, delimiter="\t" ).writerows( lstr_maf_mapping_cancer )

                        # JSON object
                        str_dna_rna_file_base = os.path.splitext( os.path.basename( str_call_run_conf ) )[ 0 ]
                        str_dna_rna_json_file = os.path.join( args_parsed.str_file_base, str_dna_rna_file_base, str_dna_rna_file_base + ".json" )
                        str_dna_rna_json_command = " ".join( [ "make_inspector_json.py",
                                                               "--output_file", str_dna_rna_json_file ] +
                                                             [ "--input_files " + str_dna_rna_file_info 
                                                               for str_dna_rna_file_info in lstr_dna_rna_file_info ] )
                        cmd_dna_rna_json = Command.Command( str_cur_command = str_dna_rna_json_command,
                                                    lstr_cur_dependencies = lstr_dna_rna_json_file_dependencies,
                                                    lstr_cur_products = str_dna_rna_json_file )
                        lcmd_commands_run.append( cmd_dna_rna_json )

                        # Figures
                        dict_validation_commands = self.func_validation_figure_commands( args_parsed = args_parsed,
                                                                        str_cur_project_dir = str_current_project_dir,
                                                                        lstr_dna_vcfs_snps = lstr_generated_dna_vcf_snps,
                                                                        lstr_rna_vcfs_snps = lstr_generated_rna_vcf_snps,
                                                                        lstr_maf_rna_tab = [],
                                                                        lstr_maf_dna_tab = [],
                                                                        lstr_dna_rna_tab = lstr_generated_dna_rna_tab,
                                                                        lstr_dna_rna_tab_cancer = lstr_generated_dna_rna_tab_cancer,
                                                                        lstr_dna_rna_tab_cosmic = lstr_generated_dna_rna_tab_cosmic,
                                                                        lstr_dna_rna_tab_edit = lstr_generated_dna_rna_tab_edit,
                                                                        lstr_syn_rna_tab = [],
                                                                        lstr_eval_init_vcf = lstr_dna_rna_dependencies,
                                                                        lstr_eval_edit_vcf = lstr_dna_rna_edit_dependencies,
                                                                        lstr_eval_cosmic_vcf = lstr_dna_rna_cosmic_dependencies, 
                                                                        lstr_eval_cancer_vcf = lstr_dna_rna_cravat_dependencies,
                                                                        str_mapping_init = str_maf_mapping_init_file,
                                                                        str_mapping_edit = str_maf_mapping_edit_file,
                                                                        str_mapping_cosmic = str_maf_mapping_cosmic_file,
                                                                        str_mapping_cancer = str_maf_mapping_cancer_file )
                        lcmd_commands_run.extend( dict_validation_commands[ STR_CMDS ] )
                    """
                    ######################################################
                    # Make synthetic comparisons tab files
                    # Make tabs for synthetic truth data and RNASEQ Data
                    #####################################################
                    str_synthetic_tab_sample_file = os.path.join( str_current_project_dir, "tab_syn_samples.txt" )
                    str_synthetic_tab_run_conf_file = os.path.join( args_parsed.str_tab_config_dir, "run_synthetic_tabs.conf" )
                    lstr_syn_tab_products = []
                    lstr_syn_link_products = []
                    llstr_synthetic_samples = []
                    lstr_syn_symbolic_link_command = []
                    lstr_syn_file_info = []
                    lstr_syn_json_file_dependencies = []
                    for dict_sample in dict_sample_study.values():
                        if not dict_sample[ STR_SYNTHETIC_VCF ]:
                            continue

                        str_cur_tab_sample = os.path.splitext( os.path.basename( dict_sample[ STR_RNA_LEFT ] ) )[0]
                        str_link_base = os.path.join( str_current_project_dir, "samples", str_cur_tab_sample, "reads." + str_cur_tab_sample )
                        str_cur_tab_file = self.func_convert_fastq_syn_rna_tab( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                        llstr_synthetic_samples.append( "\t".join( [ str_cur_tab_sample,
                                                                 os.path.join( str_current_project_dir, str_cur_tab_sample ), 
                                                                 os.path.join( str_current_project_dir, str_cur_tab_sample ) + "\n" ] ) )
                        # Link VCF files
                        str_syn_link_dep = dict_sample[ STR_SYNTHETIC_VCF ]
                        str_syn_link_prod_ref_vcf = str_link_base + "_reference_syn.vcf"
                        str_syn_link_cmd = "ln -sf " + str_syn_link_dep + " " + str_syn_link_prod_ref_vcf
                        lstr_syn_link_products.append( str_syn_link_prod_ref_vcf )
                        lcmd_commands_run.append( Command.Command( str_cur_command = str_syn_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep ],
                                                     lstr_cur_products = [ str_syn_link_prod_ref_vcf ] ) )
                        str_syn_link_prod_ref_tbi = str_link_base + "_reference_syn.vcf.tbi"
                        str_syn_index_link_cmd = "ln -sf " + str_syn_link_dep + ".tbi" + " " + str_syn_link_prod_ref_tbi
                        lcmd_commands_run.append( Command.Command( str_cur_command = str_syn_index_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep + ".tbi" ],
                                                     lstr_cur_products = [ str_syn_link_prod_ref_tbi ] ) )

                        str_syn_link_dep = self.func_convert_fastq_left_rna_vcf( dict_sample[ STR_RNA_LEFT ],str_current_project_dir )
                        str_syn_link_prod_vcf = str_link_base + "_syn.vcf"
                        str_syn_link_cmd = "ln -sf " + str_syn_link_dep + " " + str_syn_link_prod_vcf
                        lstr_syn_link_products.append( str_syn_link_prod_vcf )
                        lcmd_commands_run.append( Command.Command( str_cur_command = str_syn_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep ],
                                                     lstr_cur_products = [ str_syn_link_prod_vcf ] ) )
                        str_syn_link_prod_tbi = str_link_base + "_syn.vcf.tbi"
                        str_syn_index_link_cmd = "ln -sf " + str_syn_link_dep + ".tbi " + str_syn_link_prod_tbi
                        lcmd_commands_run.append( Command.Command( str_cur_command = str_syn_index_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep + ".tbi" ],
                                                     lstr_cur_products = [ str_syn_link_prod_tbi ] ) )

                        # Link depth files (samtools)
                        if "samtools" in str_call_run_conf.lower():
                            str_syn_link_dep = self.func_convert_fastq_rna_depth_samtools( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base )
                            str_syn_link_prod_ref_depth = str_link_base + "_reference_syn.depth"
                            str_syn_link_cmd = "ln -sf " + str_syn_link_dep + " " + str_syn_link_prod_ref_depth
                            lstr_syn_link_products.append( str_syn_link_prod_ref_depth )
                            cmd_syn_links = Command.Command( str_cur_command = str_syn_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep ],
                                                     lstr_cur_products = [ str_syn_link_prod_ref_depth ])
                            lcmd_commands_run.append( cmd_syn_links )

                            str_syn_link_dep = self.func_convert_fastq_rna_depth_samtools( dict_sample[ STR_RNA_LEFT ], args_parsed.str_file_base )
                            str_syn_link_prod_depth = str_link_base + "_syn.depth"
                            str_syn_link_cmd = "ln -sf " + str_syn_link_dep + " " + str_syn_link_prod_depth
                            lstr_syn_link_products.append( str_syn_link_prod_depth )
                            cmd_syn_links = Command.Command( str_cur_command = str_syn_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep ],
                                                     lstr_cur_products = [ str_syn_link_prod_depth ])
                            lcmd_commands_run.append( cmd_syn_links )

                        # Link depth files (GATK)
                        if "gatk" in str_call_run_conf.lower():
                            str_syn_link_dep = self.func_convert_fastq_rna_depth_gatk( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                            str_syn_link_prod_ref_depth = str_link_base + "_reference_syn.depth"
                            str_syn_link_cmd = "ln -sf " + str_syn_link_dep + " " + str_syn_link_prod_ref_depth
                            lstr_syn_link_products.append( str_syn_link_prod_ref_depth )
                            cmd_syn_links = Command.Command( str_cur_command = str_syn_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep ],
                                                     lstr_cur_products = [ str_syn_link_prod_ref_depth ])
                            lcmd_commands_run.append( cmd_syn_links )

                            str_syn_link_dep = self.func_convert_fastq_rna_depth_gatk( dict_sample[ STR_RNA_LEFT ], str_current_project_dir )
                            str_syn_link_prod_depth = str_link_base + "_syn.depth"
                            str_syn_link_cmd = "ln -sf " + str_syn_link_dep + " " + str_syn_link_prod_depth
                            lstr_syn_link_products.append( str_syn_link_prod_depth )
                            cmd_syn_links = Command.Command( str_cur_command = str_syn_link_cmd,
                                                     lstr_cur_dependencies = [ str_syn_link_dep ],
                                                     lstr_cur_products = [ str_syn_link_prod_depth ])
                            lcmd_commands_run.append( cmd_syn_links )

                        # Tab dependencies
                        lstr_syn_tab_products.append( str_cur_tab_file )
                        lstr_generated_syn_rna_tab.append( str_cur_tab_file )

                        # Add file info for the inspector's json object
                        # test_sample,test_rna.bam,test_dna.bam,test_rna.vcf,test_dna.vcf,test_comparison.tab
                        lstr_syn_file_info.append( ",".join( [ str_cur_tab_sample, 
                                                           self.func_convert_fastq_left_rna_bam( str_cur_tab_sample, args_parsed.str_file_base ),
                                                           self.func_convert_fastq_left_dna_bam( str_cur_tab_sample, args_parsed.str_file_base ),
                                                           str_syn_link_prod_vcf, str_syn_link_prod_ref_vcf, str_cur_tab_file ] ) )
                        lstr_syn_json_file_dependencies.extend( [ self.func_convert_fastq_left_rna_bam( str_cur_tab_sample, args_parsed.str_file_base ),
                                                              self.func_convert_fastq_left_dna_bam( str_cur_tab_sample, args_parsed.str_file_base ),
                                                              str_syn_link_prod_vcf, str_syn_link_prod_ref_vcf, str_cur_tab_file ] )

                    # If there synthetic samples to run make tab files.
                    if len( llstr_synthetic_samples ) > 0:
                        if not args_parsed.f_Test:
                            with open( str_synthetic_tab_sample_file, "w" ) as hndl_syn_samples:
                                hndl_syn_samples.write( "\n".join( llstr_synthetic_samples ) )
                        str_cmd_make_synthetic_tab = " ".join( [ "run_RNASEQ_pipeline_many_samples.pl --annot_conf", args_parsed.str_annot_config,
                                                "--run_conf", str_synthetic_tab_run_conf_file, "--reads_list_file", str_synthetic_tab_sample_file,
                                                "--project_base_dir", str_current_project_dir, "--memory 20 --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                                                #"--project_base_dir", str_current_project_dir, "--memory 20 --run_on_grid --num_parallel_procs", str( args_parsed.i_jobs ) ] )
                        ############### Make tab files
                        cmd_syn_tab = Command.Command( str_cur_command = str_cmd_make_synthetic_tab,
                                                       lstr_cur_dependencies = lstr_syn_link_products + [ args_parsed.str_annot_config, str_synthetic_tab_run_conf_file ],
                                                       lstr_cur_products = lstr_syn_tab_products )
                        lcmd_commands_run.append( cmd_syn_tab )

                    ############################################
                    # Make json object for inspecting the comparison
                    # Make figures for validation with a study
                    ############################################
                    if len( lstr_generated_syn_rna_tab ) > 0:

                        # JSON object
                        str_syn_file_base = os.path.splitext( os.path.basename( str_call_run_conf ) )[ 0 ]
                        str_syn_json_file = os.path.join( args_parsed.str_file_base, str_syn_file_base, str_syn_file_base + ".json" )
                        str_syn_json_command = "make_inspector_json.py --output_file " + str_syn_json_file + " " + " ".join( [ "--input_files " + str_syn_file_info for str_syn_file_info in lstr_syn_file_info ] )
                        cmd_syn_json = Command.Command( str_cur_command = str_syn_json_command,
                                                    lstr_cur_dependencies = lstr_syn_json_file_dependencies,
                                                    lstr_cur_products = str_syn_json_file )
                        lcmd_commands_run.append( cmd_syn_json )

                        # Validation figures
                        dict_validation_commands = self.func_validation_figure_commands( args_parsed = args_parsed, str_cur_project_dir = str_current_project_dir,
                                                                        lstr_dna_vcfs_snps = [],
                                                                        lstr_rna_vcfs_snps = [],
                                                                        lstr_maf_rna_tab = [],
                                                                        lstr_maf_dna_tab = [],
                                                                        lstr_dna_rna_tab = [],
                                                                        lstr_syn_rna_tab = lstr_generated_syn_rna_tab ) 

                        lcmd_commands_run.extend( dict_validation_commands[ STR_CMDS ] )
                    """
        return( lcmd_commands_run )

    def func_validation_figure_commands( self, args_parsed,
                                         str_cur_project_dir,
                                         lstr_dna_vcfs_snps = [],
                                         lstr_rna_vcfs_snps = [],
                                         lstr_maf_rna_tab = [],
                                         lstr_maf_dna_tab = [],
                                         lstr_dna_rna_tab = [],
                                         lstr_dna_rna_tab_cancer = [],
                                         lstr_dna_rna_tab_cosmic = [],
                                         lstr_dna_rna_tab_edit = [],
                                         lstr_syn_rna_tab = [],
                                         lstr_eval_init_vcf = [],
                                         lstr_eval_edit_vcf = [],
                                         lstr_eval_cosmic_vcf = [],
                                         lstr_eval_cancer_vcf = [],
                                         str_mapping_init = None,
                                         str_mapping_edit = None,
                                         str_mapping_cosmic = None,
                                         str_mapping_cancer = None ):
        """

        * return : List of commands to be ran
                 : list of command objects
        """

        # Collects commands
        lcmd_commands = []

        # Distance metric functions
        str_jaccard_distance = os.path.join( "jaccard_distance.R" )
        str_jaccard_distance_restricted = os.path.join( "jaccard_distance_restricted.R" )

        # Figure Directories
        STR_FIGURE_DIR = os.path.join( str_cur_project_dir, self.str_figure_dir )
        STR_DNA_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_dna )
        STR_RNA_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_rna )
        STR_MAF_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_maf )
        STR_SYN_FIGURE = os.path.join( str_cur_project_dir, self.str_figure_syn )

#        if( len( lstr_dna_vcfs_snps ) > 1 and len( lstr_rna_vcfs_snps ) > 1 ):
#            # Create genotype matrix and list
#            ## DNA and RNA
#            str_genotype_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype_matrix.txt" )
#            print( str_genotype_matrix )
#            str_cmd_create_genotype_matrix = "vcfs_to_genotype_matrix.py" + " --matrix " + str_genotype_matrix + " " + " ".join( lstr_rna_vcfs_snps + lstr_dna_vcfs_snps )
#            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_create_genotype_matrix,
#                                               lstr_cur_dependencies = lstr_rna_vcfs_snps + lstr_dna_vcfs_snps,
#                                               lstr_cur_products = [ str_genotype_matrix ] ) )
#            # Visualize genotype matrices
#            ## DNA and RNA, restricted by type
#            str_genotype_matrix_restricted_pdf = os.path.join( STR_FIGURE_DIR, "sample_dna_rna_genotype_matrix_restricted.pdf" )
#            str_genotype_distance_restricted_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype_restricted.dist" )
#            print( str_genotype_matrix_restricted_pdf )
#            str_cmd_genotype_matrix_dna_rna_restricted = " ".join( [ "make_dendrogram_generic.R", "--input_matrix", str_genotype_matrix, 
#                                                      "--distance_function", str_jaccard_distance_restricted,
#                                                      "--output_pdf", str_genotype_matrix_restricted_pdf,
#                                                      "--output_distance_matrix", str_genotype_distance_restricted_matrix ] ) 
#            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_dna_rna_restricted,
#                                               lstr_cur_dependencies = [ str_genotype_matrix, str_jaccard_distance_restricted ],
#                                               lstr_cur_products = [ str_genotype_matrix_restricted_pdf, str_genotype_distance_restricted_matrix ] ) )
#            ## DNA and RNA not restricted by type
#            str_genotype_matrix_pdf = os.path.join( STR_FIGURE_DIR, "sample_dna_rna_genotype_matrix.pdf" )
#            str_genotype_distance_matrix = os.path.join( args_parsed.str_file_base, "sample_dna_rna_genotype.dist" )
#            print( str_genotype_matrix_pdf )
#            str_cmd_genotype_matrix_dna_rna = " ".join( [ "make_dendrogram_generic.R", "--input_matrix", str_genotype_matrix, 
#                                                      "--distance_function", str_jaccard_distance,
#                                                      "--output_pdf", str_genotype_matrix_pdf,
#                                                      "--output_distance_matrix", str_genotype_distance_matrix ] ) 
#            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_dna_rna,
#                                               lstr_cur_dependencies = [ str_genotype_matrix, str_jaccard_distance ],
#                                               lstr_cur_products = [ str_genotype_matrix_pdf, str_genotype_distance_matrix ] ) )
#
#        if( len( lstr_rna_vcfs_snps ) > 1 ):
#            ## RNA
#            str_genotype_rna_matrix = os.path.join( args_parsed.str_file_base, "sample_rna_genotype_matrix.txt" )
#            str_cmd_create_genotype_rna_matrix = "vcfs_to_genotype_matrix.py" + " --matrix " + str_genotype_rna_matrix + " " + " ".join( lstr_rna_vcfs_snps )
#            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_create_genotype_rna_matrix,
#                                               lstr_cur_dependencies = lstr_rna_vcfs_snps,
#                                               lstr_cur_products = [ str_genotype_rna_matrix ] ) )
#        
#            # Visualize genotype matrices
#            ## RNA
#            str_genotype_rna_matrix_pdf = os.path.join( STR_FIGURE_DIR, "sample_rna_genotype_matrix.pdf" )
#            str_genotype_rna_distance_matrix = os.path.join( args_parsed.str_file_base, "sample_rna_genotype.dist" )
#            str_cmd_genotype_matrix_rna = " ".join( [ os.path.join( "make_dendrogram_generic.R" ), "--input_matrix", str_genotype_rna_matrix, 
#                                                      "--distance_function", str_jaccard_distance,
#                                                      "--output_pdf", str_genotype_rna_matrix_pdf,
#                                                      "--output_distance_matrix", str_genotype_rna_distance_matrix ] )
#            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_genotype_matrix_rna,
#                                               lstr_cur_dependencies = [ str_genotype_rna_matrix, str_jaccard_distance ],
#                                               lstr_cur_products = [ str_genotype_rna_matrix_pdf, str_genotype_rna_distance_matrix ] ) )

        # Explore false positive and negative rates
        # DNA_RnA for initial, cosmic, RNAediting, and cancer filtering
        for str_keyword, str_dir, lstr_tabs, lstr_vcfs, str_mapping_file in zip( [ "DNA_RNA_INIT", "DNA_RNA_EDIT", "DNA_RNA_COSMIC", "DNA_RNA_CANCER" ],
                                                    [ "dna_rna", "dna_rna_edit", "dna_rna_cosmic", "dna_rna_cancer" ],
                                                    [ lstr_dna_rna_tab, lstr_dna_rna_tab_edit, lstr_dna_rna_tab_cosmic, lstr_dna_rna_tab_cancer ],
                                                    [ lstr_eval_init_vcf, lstr_eval_edit_vcf, lstr_eval_cosmic_vcf, lstr_eval_cancer_vcf ],
                                                    [ str_mapping_init, str_mapping_edit, str_mapping_cosmic, str_mapping_cancer ] ):
            if len( lstr_tabs ):
                str_DNA_RNA_figures = os.path.join( STR_FIGURE_DIR, str_dir )
                for str_tab in lstr_tabs:
                    str_cur_tab_sample = str_tab.split( os.path.sep )[ -3 ]
                    str_roc_1 = os.path.join( str_DNA_RNA_figures, "_".join( [ str_cur_tab_sample, str_keyword, "roc_truth_10_pred_vary.pdf" ] ) )
                    str_roc_2 = os.path.join( str_DNA_RNA_figures, "_".join( [ str_cur_tab_sample, str_keyword, "roc_truth_vary_pred_1.pdf" ] ) )
                    str_compare_DNA_RNA_cmd = " ".join( [ "visualize_mutation_depth_tab_files.R",
                                                  "-o", str_DNA_RNA_figures,
                                                  "-k "+str_keyword,
                                                  "-s "+str_cur_tab_sample,
                                                  "--method "+str_keyword,
                                                  str_tab ] )
                    lcmd_commands.append( Command.Command( str_cur_command = str_compare_DNA_RNA_cmd,
                                                           lstr_cur_dependencies = str_tab,
                                                           lstr_cur_products = [ str_roc_1, str_roc_2 ] ) )

            # Check the percent mutations
            if args_parsed.str_maf_file and args_parsed.str_key_mutations:
                if str_mapping_file and args_parsed.str_maf_file and args_parsed.str_key_mutations:
                    str_output_compare_maf_file = os.path.join( STR_MAF_FIGURE, os.path.basename( os.path.splitext( str_mapping_file )[0]+"_confirm_mutations.pdf" ) )
                    str_compare_maf_key_cmd = " ".join([ "confirm_maf_mutations.py",
                                                     "--maf", args_parsed.str_maf_file,
                                                     "--sample", str_mapping_file,
                                                     "--key_genes",args_parsed.str_key_mutations,
                                                     str_output_compare_maf_file ])
                    lcmd_commands.append( Command.Command( str_cur_command = str_compare_maf_key_cmd,
                                                           lstr_cur_dependencies = [ args_parsed.str_maf_file, str_mapping_file ] + lstr_vcfs,
                                                           lstr_cur_products = [ str_output_compare_maf_file ]) )

            # Make vchk files for the vcf files
            if lstr_vcfs:
                lstr_vcf_cur = []
                for str_cur_vcf in  lstr_vcfs:
                    str_summary = os.path.join( STR_RNA_FIGURE, str_dir, os.path.basename( str_cur_vcf ) + "_plot.vchk" )
                    str_summary_dir = str_summary + "_dir" + os.path.sep
                    str_cur_vchk_command = " ".join([ "bcftools", "stats", str_cur_vcf, ">", str_summary ])
                    str_cur_vchk_plot_command = " ".join([ "plot-vcfstats", str_summary, "-p", str_summary_dir ])
                    lcmd_commands.append( Command.Command( str_cur_command = str_cur_vchk_command,
                                                           lstr_cur_dependencies = [ str_cur_vcf ],
                                                           lstr_cur_products = [ str_summary ] ) )
                    lcmd_commands.append( Command.Command( str_cur_command = str_cur_vchk_plot_command,
                                                           lstr_cur_dependencies = [ str_summary ],
                                                           lstr_cur_products = [ str_summary_dir ] ) )
                str_input_dir = os.path.join( STR_RNA_FIGURE, str_dir )
                str_output_dir = os.path.join( STR_RNA_FIGURE, str_dir, "group_vchk" )
                str_cmd_group_vchk = " ".join([ "combine_vchk.py", "--input_dir", str_input_dir, "--output_dir", str_output_dir ])
                lcmd_commands.append( Command.Command( str_cur_command = str_cmd_group_vchk,
                                                       lstr_cur_dependencies = [ str_input_dir ],
                                                       lstr_cur_products = [ str_output_dir ] ) )

        """ # SYNTHETIC
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

        """
        return { STR_CMDS : lcmd_commands }

if __name__ == "__main__":

    # Run pipeline
    RNASEQ_mutation_validation().func_run_pipeline()
