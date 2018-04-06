#!/usr/bin/env python


import argparse
import datetime
import os,sys
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.abspath(__file__)), "SciEDPipeR"]))
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.abspath(__file__)), "SciEDPipeR", "sciedpiper"]))
import sciedpiper.Command as Command
import sciedpiper.PipelineRunner as PipelineRunner
import sciedpiper.Pipeline as Pipeline


__author__="Timothy Tickle"
__copyright__="Copyright 2014"
__credits__=["Timothy Tickle", "Brian Haas"]
__license__="MIT"
__maintainer__="Timothy Tickle"
__email__="ttickle@broadinstitute.org"
__status__="Development"


# Constants
# Keys for alignment return (dicts)
INDEX_CMD = "cmd"
INDEX_FILE = "out_file"
INDEX_INDEX = "out_index"
INDEX_FOLDER = "align_folder"
INDEX_BAM = "bam"

# Choices for platform
STR_ILLUMINA = "ILLUMINA"
LSTR_SEQ_CHOICES = [STR_ILLUMINA, "SLX,SOLEXA",
                    "SOLID,454", "COMPLETE",
                    "PACBIO", "IONTORRENT",
                    "CAPILLARY", "HELICOS"]

# Choices for alignment
#STR_ALIGN_GSNP = "GSNAP"
STR_ALIGN_STAR = "STAR"
STR_ALIGN_STAR_LIMITED = "LIMITED"
LSTR_ALIGN_CHOICES = [STR_ALIGN_STAR, STR_ALIGN_STAR_LIMITED]

# Choices for variant calling
STR_VARIANT_GATK = "GATK"
STR_VARIANT_SAMTOOLS = "SAM"
STR_VARIANT_NONE = "NONE"
LSTR_VARIANT_CALLING_CHOICES = [STR_VARIANT_GATK,
                                STR_VARIANT_NONE]

# Choices for variant filtering
STR_FILTERING_BCFTOOLS = "BCFTOOLS"
STR_FILTERING_GATK = "GATK"
STR_FILTERING_NONE = "NONE"
STR_FILTERING_DEFAULT = STR_FILTERING_GATK
LSTR_VARIANT_FILTERING_CHOICES = [STR_FILTERING_GATK,
                                  STR_FILTERING_NONE]

# This mode is used in validating the method in the context fo DNA-seq
# It is not intended to be ran on biological samples for studies.
STR_DNASEQ_VALIDATION = "DNASEQ"

# Named files for pipeline
C_STR_CANCER_TAB = "cancer.tab"
C_STR_CANCER_ANNOTATED_VCF = "variants_annotated.vcf.gz"
C_STR_CANCER_VCF = "cancer.vcf"
C_STR_INIT_FILTER = "variants_initial_filtering.vcf"

# CRAVAT related
I_CRAVAT_ATTEMPTS = 180
I_CRAVAT_WAIT = 60
STR_CRAVAT_CLASSIFIER_DEFAULT = "Other"
STR_FDR_CUTTOFF = "0.3"

# Directory structure
STR_MISC_DIR = "misc"

# Mutation Inspector file
C_STR_MUTATION_INSPECTOR = "mutation_inspector.json"

class RnaseqSnp(PipelineRunner.PipelineRunner):
    """
    SNP detection using RNA-Seq data.
    """

    def func_do_star_alignment(self,
                               args_call,
                               str_unique_id,
                               f_index_only=False):
        """
        Manages commands for star alignment.

        * args_call : Arguments
                    : Arguments used to run pipeline
        * str_unique_id : Unique key string for this sample run,
                          to keep files unique
                        : String
        * f_index_only : Indicates that indexing only should result.
                       : Boolean (True indicates indexing only, no alignment)
        * return : Dict containing commands, resulting alignemnt file and the
                   alignment folder
               : dict
        """
        # STAR key words and staticly named files
        STR_STAR_GENOME_GENERATE = "genomeGenerate"
        STR_INDEX = "_".join(["star_index", str_unique_id])
        STR_STAR_JUNCTION = "SJ.out.tab"
        STR_STAR_LOG = "Log.out"
        STR_STAR_LOG_FINAL = "Log.final.out"
        STR_STAR_OUTPUT_BAM = "Aligned.sortedByCoord.out.bam"
        STR_STAR_PROGRESS = "Log.progress.out"

        # Pipeline variables
        str_num_threads = str(args_call.i_number_threads)
        str_left = args_call.str_sample_file_left_fq
        str_right = args_call.str_sample_file_right_fq

        # Files and dirs
        str_misc_dir = os.path.join(args_call.str_out_dir, STR_MISC_DIR)
        str_star_sorted_bam = os.path.join(str_misc_dir, STR_STAR_OUTPUT_BAM)
        str_mapped_bam_base, str_mapped_bam_ext = os.path.splitext(str_star_sorted_bam)
        str_mapped_bam = "".join([str_mapped_bam_base, "_maponly", str_mapped_bam_ext])
        str_star_output_bai = str_star_sorted_bam + ".bai"
        str_output_log_final = os.path.join(str_misc_dir,
                                            STR_STAR_LOG_FINAL)
        str_output_log = os.path.join(str_misc_dir,
                                      STR_STAR_LOG)
        str_output_log_progress = os.path.join(str_misc_dir,
                                               STR_STAR_PROGRESS)
        str_output_SJ = os.path.join(str_misc_dir,
                                     STR_STAR_JUNCTION)

        ## The index dir can be a made index or premade index.
        str_index_dir = os.path.join(args_call.str_out_dir, STR_INDEX)
        if(hasattr(args_call, "str_initial_index")
           and args_call.str_initial_index):
            str_index_dir = args_call.str_initial_index

        # Set limited memory modes
        # Update the limitGenomeGenerateRAM if more memory is requested
        lstr_index_memory_size = []
        if hasattr(args_call, "str_star_memory_limit"):
            if((not (args_call.str_star_memory_limit is None)) and
               (args_call.str_alignment_mode == STR_ALIGN_STAR)):
                lstr_index_memory_size.extend(["--limitGenomeGenerateRAM",
                                               args_call.str_star_memory_limit])

        lstr_limited_index_mode = []
        if args_call.str_alignment_mode == STR_ALIGN_STAR_LIMITED:
         lstr_limited_index_mode.append(" ".join(["--limitGenomeGenerateRAM",
                                                  "15000000000",
                                                  "--genomeSAsparseD 2",
                                                  "--outSAMmapqUnique 60",
                                                  "--limitIObufferSize",
                                                  "150000000"]))

        lstr_limited_alignment_mode = []
        if args_call.str_alignment_mode == STR_ALIGN_STAR_LIMITED:
            lstr_limited_alignment_mode.append("--genomeSAsparseD 2")

        # Commands to build and return
        lcmd_commands = []

        # If the premade index is not given then generate
        if(not hasattr(args_call, "str_initial_index")
           or not args_call.str_initial_index):
            str_index = " ".join(["STAR", "--runMode", STR_STAR_GENOME_GENERATE] +
                                  lstr_limited_index_mode + lstr_index_memory_size +
                                  ["--genomeDir", str_index_dir,
                                    "--outSAMmapqUnique", "60",
                                   "--genomeFastaFiles", args_call.str_genome_fa,
                                   "--runThreadN", str_num_threads])
            cmd_i = Command.Command(str_cur_command=str_index,
                                    lstr_cur_dependencies=[args_call.str_genome_fa],
                                    lstr_cur_products=[str_index_dir])
            lcmd_commands.append(cmd_i)

        # Perform two-pass alignment.
        if not f_index_only:

            # Star Aligner
            lstr_gzip = []
            str_ext_left = os.path.splitext(str_left)[1]
            str_ext_right = os.path.splitext(str_right)[1]
            if str_ext_left != str_ext_right:
                str_error = " ".join(["Fastq files from a single sample should both",
                                      "either be gzipped (.gz) or not compressed."])
                cur_pipeline.logr_logger.error(str_error)
                exit(8)
            elif str_ext_left == ".gz":
                lstr_gzip = ["--readFilesCommand", "\"gunzip -c\""]

            # Map files
            str_star_align = " ".join(["STAR",
                                       "--genomeDir", str_index_dir,
                                       "--runThreadN", str_num_threads,
                                       "--readFilesIn", str_left,
                                       str_right] + lstr_gzip + ["--outSAMtype",
                                       "BAM", "SortedByCoordinate",
                                       "--twopassMode", "Basic",
                                       "--limitBAMsortRAM", "30000000000",
                                       " --outSAMmapqUnique","60",
                                       "--outFileNamePrefix",
                                       str_misc_dir + os.sep])
            lstr_deps = [str_index_dir, str_left, str_right]
            lstr_prods = [str_star_sorted_bam, str_output_log_final,
                          str_output_log, str_output_log_progress, str_output_SJ]
            cmd_star_align = Command.Command(str_cur_command=str_star_align,
                                             lstr_cur_dependencies=lstr_deps,
                                             lstr_cur_products=lstr_prods)
            cmd_star_align.func_set_dependency_clean_level([str_index_dir,
                                                            str_left, str_right],
                                                            Command.CLEAN_NEVER)
            lcmd_commands.append(cmd_star_align)

            # Create bai
            str_star_bai = " ".join(["samtools index", str_star_sorted_bam])
            cmd_bai = Command.Command(str_cur_command=str_star_bai,
                                      lstr_cur_dependencies=[str_star_sorted_bam],
                                      lstr_cur_products=[str_star_output_bai])
            cmd_bai.func_set_dependency_clean_level([str_star_sorted_bam],
                                                     Command.CLEAN_NEVER)
            lcmd_commands.append(cmd_bai)

        return({INDEX_CMD:lcmd_commands,
                INDEX_FILE:str_star_sorted_bam,
                INDEX_FOLDER:str_misc_dir})

    def func_do_BWA_alignment(self,
                              args_call,
                              str_unique_id,
                              f_index_only = False):
        """
        Manages the calls for aligning DNA-seq data with BWA.
        Development use only, not intended to be used with studies.
        Used to validate technical properties of the RNA-seq pipeline.
        * args_call : Arguments
                    : Arguments used to run pipeline
        * str_unique_id : Unique key string for this sample run, to keep files unique
                        : String
        * f_index_only : Indicates that indexing only should result.
                       : Boolean (True indicates indexing only, no alignment)
        * return : Dict containing commands, resulting alignemnt file and the alignment folder
               : dict
        """
        # Files
        str_left_file_key = os.path.basename(os.path.splitext(args.str_sample_file_left_fq)[0])
        str_sam = os.path.join(args_call.str_out_dir, ".".join([str_left_file_key, "sam"]))
        str_bam = os.path.join(args_call.str_out_dir, ".".join([str_left_file_key, "bam"]))
        str_bwa_sorted_bam = os.path.join(args_call.str_out_dir, ".".join([str_left_file_key, "sorted", "bam"]))
        str_temp_prefix = os.path.join(args_call.str_out_dir, "temp")
        str_bai = str_bwa_sorted_bam + ".bai"

        lcmd_dna_mapping_commands = []

        # Pre-processing
        ## Mapping and Dedupping
        ### BWA, make coordinate ordered bam
        #### Index reference
        if f_index_only or not hasattr(args_call, "str_initial_index") or not args_call.str_initial_index:
            cmd_bwa_index = Command.Command(str_cur_command = "".join(["bwa index -a bwtsw ", args.str_genome_fa]),
                                                lstr_cur_dependencies = [args.str_genome_fa],
                                                lstr_cur_products = [args.str_genome_fa + ".amb",
                                                                     args.str_genome_fa + ".ann",
                                                                     args.str_genome_fa + ".bwt",
                                                                     args.str_genome_fa + ".pac",
                                                                     args.str_genome_fa + ".sa"])
            if f_index_only:
                return { INDEX_CMD: [cmd_bwa_index], INDEX_FOLDER: os.path.dirname(args.str_genome_fa) }
            else:
                lcmd_dna_mapping_commands.append(cmd_bwa_index)

        # Align both samples
        # bwa sample ref.fasta fwd.sai rev.sai fwd.fq rev.fq > mydata.sam
        cmd_bwa_sam = Command.Command(str_cur_command = "".join(["bwa mem -M -R \"@RG\\tID:", str_unique_id,"\\tSM:", str_unique_id,"\" ",
                                                                   args.str_genome_fa, " ", args.str_sample_file_left_fq, " ", args.str_sample_file_right_fq, " > ", str_sam]),
                                                lstr_cur_dependencies = [args.str_genome_fa, args.str_sample_file_left_fq, args.str_sample_file_right_fq],
                                                lstr_cur_products = [str_sam])
        lcmd_dna_mapping_commands.append(cmd_bwa_sam)

        # SAM to BAM
        lcmd_dna_mapping_commands.append(Command.Command(str_cur_command = " ".join(["samtools view -b -S -o", str_bam, str_sam]),
                                               lstr_cur_dependencies = [str_sam],
                                               lstr_cur_products = [str_bam]))

        # Sort coordinate order
        lcmd_dna_mapping_commands.append(Command.Command(str_cur_command = " ".join(["samtools sort -O bam -T " + str_temp_prefix + " -o", str_bwa_sorted_bam, str_bam]),
                                               lstr_cur_dependencies = [str_bam],
                                               lstr_cur_products = [str_bwa_sorted_bam]))

        # Create bai
        lcmd_dna_mapping_commands.append(Command.Command(str_cur_command = " ".join(["samtools index", str_bwa_sorted_bam]),
                                               lstr_cur_dependencies = [str_bwa_sorted_bam],
                                               lstr_cur_products = [str_bai]))
        return { INDEX_CMD: lcmd_dna_mapping_commands, INDEX_FILE: str_bwa_sorted_bam, INDEX_FOLDER:args_call.str_out_dir }

    def func_do_recalibration_gatk(self, args_call, str_align_file,
                                   str_unique_id, str_project_dir, str_tmp_dir,
                                   lstr_dependencies, logr_cur):
        """
        Manages the commands for the recalibration step in the GATK RNASEq mutation calling best practices.

        * args_call : Arguments for the pipeline
                    : Dict
        * str_align_file : The file from the alignment (sam or bam file). If sam file, will be changed to bam file
                         : File path
        * str_unique_id : Key string for the smaple run (to keep files unique)
                        : String
        * str_project_dir : Output directory
                          : String file path
        * str_tmp_dir : Directory used to put intermediary files (part of the pipeline organization
                      : String file path
        * lstr_dependencies : List of file paths of dependencies from any previously running commands.
                            : List of strings
        * logr_cur : Pipeline logger
                   : Logger
        """
        # Check for the known vcf file
        # If it does not exist, warn that the associated steps will not be ran.
        if args_call.str_vcf_file is None:
            logr_cur.warn("".join(["\n\n\nWARNING, WARNING, WARNING, WARNING.\n",
                          "GATK Recalibration: A vcf file with known variants was not provided for realignment and recalibration steps.\n",
                           "These steps may perform better if such a vcf file is provided.\n\n\n"]))

        # Files
        str_dedupped_bam = os.path.join(str_tmp_dir, "dedupped.bam")
        str_dedupped_bai = os.path.join(str_tmp_dir, "dedupped.bai")
        str_intervals = os.path.join(str_tmp_dir, "forIndelRealigner.intervals")
        str_qc_metrics = os.path.join(str_tmp_dir, "mark_duplicates_qc_metrics.txt")
        str_realigned_bam = os.path.join(str_tmp_dir, "realigned.bam")
        str_realigned_bai = os.path.join(str_tmp_dir, "realigned.bai")
        str_recalibrated_alignment_file = os.path.join(str_tmp_dir, "recal_table.table")
        str_recalibrated_bam = os.path.join(str_tmp_dir, "recalibrated.bam")
        str_recalibrated_bam_2 = os.path.join(str_tmp_dir, "recalibrated_tmp.bam")
        str_recalibrated_bai = os.path.join(str_tmp_dir, "recalibrated.bai")
        str_recalibrated_bai_2 = os.path.join(str_tmp_dir, "recalibrated_tmp.bai")
        str_recalibration_plots_pdf = os.path.join(str_tmp_dir, "recalibration.pdf")
        str_sorted_bam = os.path.join(str_tmp_dir, "sorted.bam")
        str_split_bam = os.path.join(str_tmp_dir, "split.bam")
        str_split_bai = os.path.join(str_tmp_dir, "split.bai")
        # This is the file that is returned, could be many of the files below depending on the settings
        # This is dynamically set and different parts of this pipeline segment are activated.
        str_return_bam = ""
        # Allows the known variants vcf file to be available or not.
        lstr_known_vcf = [] if args_call.str_vcf_file is None else ["-known", args_call.str_vcf_file]
        lstr_known_two_dash_vcf = [] if args_call.str_vcf_file is None else ["--known", args_call.str_vcf_file]
        lcmd_gatk_recalibration_commands = []
        # Create commands
        # SAM to BAM and QC
        cmd_add_or_replace_groups = Command.Command(str_cur_command = "".join(["java -jar AddOrReplaceReadGroups.jar I=", str_align_file,
                                                                         " O=", str_sorted_bam, " SO=coordinate RGID=id RGLB=library RGPL=",
                                                                         args_call.str_sequencing_platform, " RGPU=machine RGSM=", str_unique_id]),
                                                lstr_cur_dependencies = lstr_dependencies,
                                                lstr_cur_products = [str_sorted_bam])
        cmd_add_or_replace_groups.func_set_dependency_clean_level([str_align_file], Command.CLEAN_NEVER)
        cmd_mark_duplicates = Command.Command(str_cur_command = "".join(["java -jar MarkDuplicates.jar I=", str_sorted_bam, " O=", str_dedupped_bam,
                                                                          " CREATE_INDEX=true M=", str_qc_metrics]),
                                                lstr_cur_dependencies = [str_sorted_bam],
                                                lstr_cur_products = [str_dedupped_bam, str_dedupped_bai, str_qc_metrics])
        cmd_split_cigar_reads = Command.Command(str_cur_command = " ".join(["java -jar gatk-package-4.0.1.2-local.jar SplitNCigarReads -R", args_call.str_genome_fa,
                                                                          "-I", str_dedupped_bam, "-O", str_split_bam,
                                                                          "--read-validation-stringency LENIENT"]),
                                                lstr_cur_dependencies = [args_call.str_genome_fa, str_dedupped_bam, str_dedupped_bai],
                                                lstr_cur_products = [str_split_bam, str_split_bai])
        lcmd_gatk_recalibration_commands.extend([cmd_add_or_replace_groups, cmd_mark_duplicates, cmd_split_cigar_reads])
        str_return_bam = str_split_bam
        str_return_bai = str_split_bai
        '''

        # optional indel realignment step
        if not args_call.f_stop_optional_realignment:
            cmd_realigner_target_creator = Command.Command(str_cur_command = " ".join(["java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R",
                                                                              args_call.str_genome_fa, "-I", str_split_bam, "--out", str_intervals] + lstr_known_two_dash_vcf),
                                                    lstr_cur_dependencies = [args_call.str_genome_fa, str_split_bam, str_split_bai, args_call.str_vcf_file],
                                                    lstr_cur_products = [str_intervals])
            cmd_indel_realigner = Command.Command(str_cur_command = " ".join(["java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R",
                                                                               args_call.str_genome_fa, "-I", str_split_bam, "-targetIntervals",
                                                                               str_intervals, "--out", str_realigned_bam] + lstr_known_vcf),
                                                     lstr_cur_dependencies = [args_call.str_genome_fa, str_split_bam, str_split_bai, str_intervals, args_call.str_vcf_file],
                                                     lstr_cur_products = [str_realigned_bam, str_realigned_bai])
            lcmd_gatk_recalibration_commands.extend([cmd_realigner_target_creator, cmd_indel_realigner])
            str_return_bam = str_realigned_bam
            str_return_bai = str_realigned_bai
        '''
        # Recalibrate alignments
        if not args_call.str_vcf_file is None:
            cmd_base_recalibrator = Command.Command(str_cur_command = " ".join(["java -Xmx4g -jar gatk-package-4.0.1.2-local.jar BaseRecalibrator -I",
                                                                          str_split_bam,
                                                                          "-R", args_call.str_genome_fa, "-O", str_recalibrated_alignment_file, "--known-sites", args_call.str_vcf_file]),
                                                lstr_cur_dependencies = [args_call.str_genome_fa, args_call.str_vcf_file] +
                                                                      [str_split_bam, str_split_bai],
                                                lstr_cur_products = [str_recalibrated_alignment_file])
           
            cmd_print_reads = Command.Command(str_cur_command = " ".join(["java -Xmx2g -jar gatk-package-4.0.1.2-local.jar ",
                                                                           "PrintReads", "-O", str_recalibrated_bam_2, "-I",
                                                                           str_split_bam]),
                                                 lstr_cur_dependencies = [args_call.str_genome_fa, str_recalibrated_alignment_file] +
                                                                      [str_split_bam, str_split_bai] ,
                                                 lstr_cur_products = [str_recalibrated_bam_2,str_recalibrated_bai_2])

            cmd_apply_bqsr = Command.Command(str_cur_command = " ".join(["java -jar gatk-package-4.0.1.2-local.jar ApplyBQSR",
                                                    "-I", str_recalibrated_bam_2, "-O", str_recalibrated_bam,
                                                    "-bqsr", str_recalibrated_alignment_file]),
                                                 lstr_cur_dependencies = [args_call.str_genome_fa, str_recalibrated_alignment_file] +
                                                                      [str_split_bam, str_split_bai] ,
                                                 lstr_cur_products = [str_recalibrated_bam, str_recalibrated_bai])
            
            lcmd_gatk_recalibration_commands.extend([cmd_base_recalibrator, cmd_print_reads, cmd_apply_bqsr])
            str_return_bam = str_recalibrated_bam
            str_return_bai = str_recalibrated_bai

            # Optional plotting of recalibration
            if args_call.f_optional_recalibration_plot:
                cmd_analyse_covariates = Command.Command(str_cur_command = " ".join(["java -Xmx4g -jar gatk-package-4.0.1.2-local.jar AnalyzeCovariates -R",
                                                                            args_call.str_genome_fa, "-bqsr", str_recalibrated_alignment_file,
                                                                            "-plots", str_recalibration_plots_pdf]),
                                                  lstr_cur_dependencies = [args_call.str_genome_fa, str_recalibrated_alignment_file],
                                                  lstr_cur_products = [str_recalibration_plots_pdf])
                cmd_analyse_covariates.func_set_dependency_clean_level([str_recalibrated_alignment_file], Command.CLEAN_NEVER)
                lcmd_gatk_recalibration_commands.append(cmd_analyse_covariates)
        return {INDEX_CMD:lcmd_gatk_recalibration_commands,
                INDEX_FILE:str_return_bam,
                INDEX_INDEX:str_return_bai,
                INDEX_FOLDER:str_tmp_dir}

    def func_do_rnaseq_caller_gatk(self, args_call,
                                   str_input_bam,
                                   str_input_bai,
                                   str_unique_id,
                                   str_project_dir,
                                   str_tmp_dir):
        """
        Creates the commands for the GATK RNASeq calling.

        * args_call : Arguments for running the pipeline
                    : dict
        * str_input_bam : bam to call against
                        : string file path
        * str_unique_id : Unique key for this sample run used to keep files unique
                        : string
        * str_project_dir : Output directory for project
                          : string path
        * str_tmp_dir : Directory used to place intermediary files
                        in the pipeline (mainly for organization)
                      : stribg directory
        * return : List of commands
                 : list
        """
        # Commands
        lcmd_gatk_rna_calling = []

        # Files
        str_variants_file = os.path.join(str_project_dir, "variants.vcf")

        # Create depth file
        if args_call.f_calculate_base_coverage:
    #        str_depth_compressed_file = os.path.basename(args_call.str_out_dir) + ".depth.gz"
            str_depth_compressed_file = os.path.basename(args_call.str_out_dir) + ".depth"
            str_depth_compressed_file = os.path.join(args_call.str_out_dir,
                                                     str_depth_compressed_file)
    #        lcmd_samtools_variants_commands.append(Command.Command(str_cur_command = "samtools depth " + str_input_bam + " | gzip > " + str_depth_compressed_file,
            str_depth = " ".join(["samtools", "depth",
                                  str_input_bam, ">", str_depth_compressed_file])
            cmd_depth = Command.Command(str_cur_command=str_depth,
                                        lstr_cur_dependencies=[str_input_bam],
                                        lstr_cur_products=[str_depth_compressed_file])
            cmd_depth.func_set_dependency_clean_level([str_input_bam], Command.CLEAN_NEVER)
            lcmd_gatk_rna_calling.append(cmd_depth)

        # Variant calling
        str_hap_call = " ".join(["java", "-jar", "gatk-package-4.0.1.2-local.jar",
                                 "HaplotypeCaller", "-R",
                                 args_call.str_genome_fa, "-I", str_input_bam,
                                 "--recover-dangling-heads","true",
                                 "--dont-use-soft-clipped-bases","-stand-call-conf",
                                 "20.0",
                                 "-O", str_variants_file])
        cmd_haplotype_caller = Command.Command(str_cur_command=str_hap_call,
                                               lstr_cur_dependencies=[args_call.str_genome_fa,
                                                                      str_input_bam,
                                                                      str_input_bai],
                                               lstr_cur_products=[str_variants_file])
        cmd_haplotype_caller.func_set_dependency_clean_level([str_input_bam,
                                                              str_input_bai],
                                                             Command.CLEAN_NEVER)
        lcmd_gatk_rna_calling.append(cmd_haplotype_caller)

        return {INDEX_CMD:lcmd_gatk_rna_calling,
                INDEX_FILE:str_variants_file}

    def func_do_variant_calling_gatk(self, args_call,
                                     str_align_file,
                                     str_unique_id,
                                     str_project_dir,
                                     str_tmp_dir,
                                     lstr_dependencies,
                                     logr_cur):
        """
        Creates the commands for the GATK variant calling pipeline.

        * args_call: Arguments for the pipeline
                   : Dict
        * str_align_file: The file from the alignment (sam or bam file).
                          If sam file, will be changed to bam file
                        : File path
        * str_unique_id: Key string for the smaple run (to keep files unique)
                       : String
        * str_project_dir: Output directory
                         : String file path
        * str_tmp_dir: Directory used to put intermediary files
                       (part of the pipeline organization
                     : String file path
        * lstr_dependencies: List of file paths of dependencies from
                             any previously running commands.
                           : List of strings
        * logr_cur: Pipeline logger
                  : Logger
        * return: List of commands
                : List of commands to run for BWA alignment
        """

        # Commands which will be returned
        lcmd_gatk_variants_commands = []

        # Perform recalibration
        dict_recal = self.func_do_recalibration_gatk(args_call,
                                                str_align_file,
                                                str_unique_id,
                                                str_project_dir,
                                                str_tmp_dir,
                                                lstr_dependencies,
                                                logr_cur)
        lcmd_gatk_variants_commands.extend(dict_recal[INDEX_CMD])

        # Do calling
        dict_rnaseq_gatk = self.func_do_rnaseq_caller_gatk(args_call,
                                                      dict_recal[INDEX_FILE],
                                                      dict_recal[INDEX_INDEX],
                                                      str_unique_id,
                                                      str_project_dir,
                                                      str_tmp_dir)
        lcmd_gatk_variants_commands.extend(dict_rnaseq_gatk[INDEX_CMD] )

        return({INDEX_CMD:lcmd_gatk_variants_commands,
                INDEX_FILE:dict_rnaseq_gatk[INDEX_FILE],
                INDEX_BAM:dict_recal[INDEX_FILE],
                INDEX_INDEX:dict_recal[INDEX_INDEX]})

    def func_call_dnaseq_like_rnaseq(self, args_call, str_align_file,
                                     str_unique_id, str_project_dir,
                                     str_tmp_dir, lstr_dependencies, logr_cur):
        """
        Manages the calls for calling mutations in DNA-seq data in a similar way to the RNA-seq calls.
        Development use only, not intended to be used with studies.
        Used to validate technical properties of the RNA-seq pipeline.

        Creates the commands for the SamTools variant calling pipeline.

        * args_call : Arguments for the pipeline
                    : Dict
        * str_align_file : The file from the alignment (sam or bam file). If sam file, will be changed to bam file
                         : File path
        * str_unique_id : Key string for the smaple run (to keep files unique)
                        : String
        * str_project_dir : Output directory
                          : String file path
        * str_tmp_dir : Directory used to put intermediary files (part of the pipeline organization
                      : String file path
        * lstr_dependencies : List of file paths of dependencies from any previously running commands.
                            : List of strings
        * logr_cur : Pipeline logger
                   : Logger
        * return : List of commands
               : List of commands to run for BWA alignment
        """
        str_sorted_bam = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.bam"]))
        str_sorted_bam_bai = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.bam.bai"]))
        str_dedup_bam = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.dedupped.bam"]))
        str_dedup_metrics = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.dedup.metrics"]))
        str_filtered_variants_file = os.path.join(str_project_dir, ".".join([str_unique_id, "filtered.variants.vcf"]))
        str_intervals = os.path.join(str_tmp_dir, ".".join([str_unique_id, "religner.intervals"]))
        str_raw_vcf = os.path.join(str_tmp_dir, ".".join([str_unique_id, "variants.vcf"]))
        str_realigned_bam = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.dedup.groups.realigned.bam"]))
        str_realigned_bai = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.dedup.groups.realigned.bai"]))
        str_recal_plot = os.path.join(str_tmp_dir, ".".join([str_unique_id, "recal.pdf"]))
        str_recal_snp_bam = os.path.join(str_tmp_dir, ".".join([str_unique_id, "recal_snp.bam"]))
        str_recal_snp_bai = os.path.join(str_tmp_dir, ".".join([str_unique_id, "recal_snp.bai"]))
        str_recal_table = os.path.join(str_tmp_dir, ".".join([str_unique_id, "recal.table"]))
        str_recal_table_2 = os.path.join(str_tmp_dir, ".".join([str_unique_id, "recal_2.table"]))
        str_replace_bam = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.dedup.groups.bam"]))
        str_replace_bai = os.path.join(str_tmp_dir, ".".join([str_unique_id, "sorted.dedup.groups.bai"]))

        # DNA-seq best practices
        # java -jar SortSam.jar I=Input.sam O=output.bam SO=coordinate
        cmd_sort_bam = Command.Command(str_cur_command = "".join(["java -jar SortSam.jar SO=coordinate I=", str_align_file, " O=", str_sorted_bam]),
                                                lstr_cur_dependencies = [str_align_file],
                                                lstr_cur_products = [str_sorted_bam])

        # Create bai
        cmd_sort_index_bam = Command.Command(str_cur_command = " ".join(["samtools index", str_sorted_bam]),
                              lstr_cur_dependencies = [str_sorted_bam],
                              lstr_cur_products = [str_sorted_bam_bai])
        cmd_sort_index_bam.func_set_dependency_clean_level([str_sorted_bam], Command.CLEAN_NEVER)

        # java -jar MarkDuplicates.jar I=input.sam O=output.bam
        cmd_dedup = Command.Command(str_cur_command = "".join(["java -jar MarkDuplicates.jar I=", str_sorted_bam,
                                                                 " M=", str_dedup_metrics, " O=", str_dedup_bam]),
                                                lstr_cur_dependencies = [str_sorted_bam],
                                                lstr_cur_products = [str_dedup_metrics, str_dedup_bam])

        # java -jar AddOrReplaceReadGroups.jar I=input.bam O=output.bam RGID=x RGLB=x RGPL=x RGPU=x RGSM=x RGCN=x RGDT=x
        cmd_replace = Command.Command(str_cur_command = "".join(["java -jar AddOrReplaceReadGroups.jar I=", str_dedup_bam, " O=", str_replace_bam,
                                                                  " RGCN=RGCN RGID=", str_unique_id, " RGLB=library RGDT=", datetime.date.today().isoformat()," RGPL=",
                                                                    args_call.str_sequencing_platform, " RGPU=machine RGSM=", str_unique_id, " CREATE_INDEX=TRUE"]),
                                                lstr_cur_dependencies = [str_dedup_bam],
                                                lstr_cur_products = [str_replace_bam, str_replace_bai])

        # Commands so far
        ls_cmds = [cmd_sort_bam, cmd_sort_index_bam, cmd_dedup, cmd_replace, cmd_create_target, cmd_realign, cmd_recalibrate, cmd_print, cmd_recalibrate_2]
        # Create depth file
        if args_call.f_calculate_base_coverage:
    #        str_depth_compressed_file = os.path.basename(args_call.str_out_dir) + ".depth.gz"
            str_depth_compressed_file = os.path.basename(args_call.str_out_dir) + ".depth"
            str_depth_compressed_file = os.path.join(args_call.str_out_dir, str_depth_compressed_file)
    #        lcmd_samtools_variants_commands.append(Command.Command(str_cur_command = "samtools depth " + str_recal_snp_bam + " | gzip > " + str_depth_compressed_file,
            ls_cmds.append(Command.Command(str_cur_command = "samtools depth " + str_recal_snp_bam + " > " + str_depth_compressed_file,
                                                   lstr_cur_dependencies = [str_recal_snp_bam],
                                                   lstr_cur_products = [str_depth_compressed_file]))

        # Call mutations - Single File, variant only calling in DNA-seq.
        # https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        # java -jar GenomeAnalysisTk.jar -T HaplotypeCaller -R reference/file.fasta -I recal.bam -stand_call_conf 30 -stand_emit_conf -o output.vcf
        cmd_haplotype_caller = Command.Command(str_cur_command = " ".join(["java -jar gatk-package-4.0.1.2-local.jar HaplotypeCaller -R", args_call.str_genome_fa,
                                                               "-I", str_recal_snp_bam, "-stand-call-conf 30.0 -O", str_raw_vcf]),
                                                lstr_cur_dependencies = [args_call.str_genome_fa, str_recal_snp_bam, str_recal_snp_bai],
                                                lstr_cur_products = [str_raw_vcf])
        cmd_haplotype_caller.func_set_dependency_clean_level([str_recal_snp_bam, str_recal_snp_bai], Command.CLEAN_NEVER)

        ls_cmds.extend([cmd_haplotype_caller]) #, cmd_variant_filteration])
        return {INDEX_CMD:ls_cmds,
                INDEX_FILE:str_raw_vcf,
                INDEX_BAM:str_recal_snp_bam,
                INDEX_INDEX:str_recal_snp_bai}

    def func_do_variant_calling_samtools(self, args_call, str_align_file,
                                         str_unique_id, str_project_dir,
                                         str_tmp_dir, lstr_dependencies,
                                         logr_cur):
        """
        Creates the commands for the SamTools variant calling pipeline.

        * args_call : Arguments for the pipeline
                    : Dict
        * str_align_file : The file from the alignment (sam or bam file). If sam file, will be changed to bam file
                         : File path
        * str_unique_id : Key string for the smaple run (to keep files unique)
                        : String
        * str_project_dir : Output directory
                          : String file path
        * str_tmp_dir : Directory used to put intermediary files (part of the pipeline organization
                      : String file path
        * lstr_dependencies : List of file paths of dependencies from any previously running commands.
                            : List of strings
        * logr_cur : Pipeline logger
                   : Logger
        """
        # Commands to run
        lcmd_samtools_variants_commands = []
        # The bam file, stored here because the path may be changed if the optional sam -> conversion is needed.
        str_bam_sorted = str_align_file
        # Index of the sorted bam file
        str_bam_sorted_index = ".".join([str_bam_sorted, "bai"])
        # Uncompressed variant calling file
        str_variants_vcf = os.path.join(str_tmp_dir, ".".join([str_unique_id, "vcf"]))

        # Optional SAM to BAM
        if os.path.splitext(str_align_file)[1].lower() == ".sam":
            str_bam = ".".join([os.path.splitext(str_align_file)[0],"bam"])
            # Sorted bam file path
            str_bam_file = os.path.split(str_bam)[1]
            str_bam_sorted = os.path.join(str_tmp_dir, self.func_switch_ext(str_bam_file, "_sorted.bam"))
            str_temp_prefix = os.path.join(str_tmp_dir, temp)
            str_bam_sorted_index = ".".join([str_bam_sorted, "bai"])
            lcmd_samtools_variants_commands.extend([
                                Command.Command(str_cur_command = " ".join(["samtools view -b -S -o",str_bam, str_align_file]),
                                                lstr_cur_dependencies = lstr_dependencies,
                                                lstr_cur_products = [str_bam]),
                                Command.Command(str_cur_command = " ".join(["samtools sort -O bam -T " + str_temp_prefix + " -o", str_bam_sorted, str_bam]),
                                                lstr_cur_dependencies = [str_bam],
                                                lstr_cur_products = [str_bam_sorted]),
                                Command.Command(str_cur_command = " ".join(["samtools index", str_bam_sorted]),
                                                lstr_cur_dependencies = [str_bam_sorted],
                                                lstr_cur_products = [str_bam_sorted_index])])

        # Either prepare bams with GATK best practices or minimally
        if args_call.f_recalibrate_sam:
            # Create commands for recalibration
            # Update files to recalibrated files
            dict_recalibration = self.func_do_recalibration_gatk(args_call, str_bam_sorted, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies, logr_cur)
            lcmd_samtools_variants_commands.extend(dict_recalibration[INDEX_CMD])
            str_bam_sorted = dict_recalibration[INDEX_FILE]
            str_bam_sorted_index = dict_recalibration[INDEX_INDEX]

        # Create depth file
        if args_call.f_calculate_base_coverage:
    #        str_depth_compressed_file = os.path.basename(args_call.str_out_dir + ".depth.gz")
            str_depth_compressed_file = os.path.basename(args_call.str_out_dir + ".depth")
            str_depth_compressed_file = os.path.join(args_call.str_out_dir, str_depth_compressed_file)
    #        lcmd_samtools_variants_commands.append(Command.Command(str_cur_command = "samtools depth " + str_bam_sorted + " | gzip > " + str_depth_compressed_file,

            cmd_depth_file = Command.Command(str_cur_command = "samtools depth " + str_bam_sorted + " > " + str_depth_compressed_file,
                                              lstr_cur_dependencies = [str_bam_sorted],
                                              lstr_cur_products = [str_depth_compressed_file])
            cmd_depth_file.func_set_dependency_clean_level([str_bam_sorted], Command.CLEAN_NEVER)
            lcmd_samtools_variants_commands.append(cmd_depth_file)

        # Identify variants
        str_samtools_calls = " ".join(["samtools mpileup -ugf", args_call.str_genome_fa, str_bam_sorted, "| bcftools call -mv -Ov >", str_variants_vcf])
        cmd_sam_call =  Command.Command(str_cur_command = str_samtools_calls,
                                         lstr_cur_dependencies = [str_bam_sorted],
                                         lstr_cur_products = [str_variants_vcf])
        cmd_sam_call.func_set_dependency_clean_level([str_bam_sorted], Command.CLEAN_NEVER)
        lcmd_samtools_variants_commands.append(cmd_sam_call)
        return {INDEX_CMD:lcmd_samtools_variants_commands,
                INDEX_FILE:str_variants_vcf,
                INDEX_BAM:str_bam_sorted,
                INDEX_INDEX:str_bam_sorted_index}

    def func_do_variant_calling_none(self, args_call, str_align_file,
                                     str_unique_id, str_project_dir,
                                     str_tmp_dir, lstr_dependencies, logr_cur):
        """
        Does nothing, allows no calling
        """
        return{INDEX_CMD:[],
               INDEX_FILE:"",
               INDEX_BAM:"",
               INDEX_INDEX:""}

    def func_do_filtering_bcftools(self, args_call, str_variants_file,
                                   lstr_dependencies, logr_cur):
        """
        Creates the commands for the bcftools hard filtering
        and custom variant cluster filtering.

        * args_call : Arguments for the pipeline
                    : Dict
        * str_variants_file : Path to file to be filtered
                            : String path
        * lstr_dependencies : List of file paths of dependencies from any
                              previously running commands.
                            : List of strings
        * logr_cur : Pipeline logger
                   : Logger
        """
        # Filtered variants file
        str_standard_variants_file = args_call.str_out_dir + os.sep + C_STR_INIT_FILTER
        str_filtered_variants_file = self.func_switch_ext(str_variants_file, "_decluster.vcf")
        str_filtered_variants_index_file = str_filtered_variants_file + ".csi"

        # Filter variants
        str_filter_command = " ".join(["bcftools", "filter", "--output-type v",
                                       "--output", str_standard_variants_file,
                                       "-sLowQual","-g 3","-G 10",
                                       "-e \"%QUAL<10 || (RPB<0.1 && %QUAL<15)",
                                       "|| (AC<2 && %QUAL<15)\"",
                                       str_variants_file])
        cmd_variant_filtration = Command.Command(str_cur_command=str_filter_command,
                                                lstr_cur_dependencies=[args_call.str_genome_fa] + lstr_dependencies,
                                                lstr_cur_products=[str_standard_variants_file])

        # Filter out clusters of SNPs
        str_custom_filter_command = " ".join(["filter_variant_clusters.py",
                                              "--window 35 --cluster 3",
                                              str_standard_variants_file,
                                              str_filtered_variants_file])
        cmd_secondary_filters = Command.Command(str_cur_command=str_custom_filter_command,
                                                lstr_cur_dependencies=[str_standard_variants_file],
                                                lstr_cur_products=[str_filtered_variants_file])

        return {INDEX_CMD:[cmd_variant_filtration, cmd_secondary_filters],
                INDEX_FILE:str_filtered_variants_file}

    def func_do_filtering_gatk(self, args_call, str_variants_file,
                               lstr_dependencies, logr_cur):
        """
        Creates the commands for the gatk hard filtering.
        * args_call: Arguments for the pipeline
                   : Dict
        * str_variants_file: Path to file to be filtered
                           : String path
        * lstr_dependencies: List of file paths of dependencies from any
                             previously running commands.
                           : List of strings
        * logr_cur: Pipeline logger
                  : Logger
        """
        # Filtered variants file
        str_filtered_variants_file = os.path.join(args_call.str_out_dir,
                                                  C_STR_INIT_FILTER)
        str_filtered_variants_index_file = str_filtered_variants_file + ".csi"
        # Filter variants
        str_filter_command = " ".join(["java -jar gatk-package-4.0.1.2-local.jar",
                                       "VariantFiltration -R",
                                       args_call.str_genome_fa, "-V",
                                       str_variants_file, "-window 35",
                                       "-cluster 3 --filter-name FS",
                                       "-filter \"FS > 30.0\"",
                                       "--filter-name QD","-filter \"QD < 2.0\"",
                                       "-O", str_filtered_variants_file])
        cmd_variant_filteration = Command.Command(str_cur_command=str_filter_command,
                                                   lstr_cur_dependencies=[args_call.str_genome_fa] + lstr_dependencies,
                                                   lstr_cur_products=[str_filtered_variants_file])

        return{INDEX_CMD: [cmd_variant_filteration],
               INDEX_FILE: str_filtered_variants_file}

    def func_do_filtering_none(self, args_call, str_variants_file,
                               lstr_dependencies, logr_cur):
        """
        Creates the commands for the gatk hard filtering.
        * args_call: Arguments for the pipeline
                   : Dict
        * str_variants_file: Path to file to be filtered
                           : String path
        * lstr_dependencies: List of file paths of dependencies from
                             any previously running commands.
                           : List of strings
        * logr_cur: Pipeline logger
                  : Logger
        """
        return { INDEX_CMD : [], INDEX_FILE : "" }

    def func_do_variant_filtering_cancer(self, args_call, str_variants_file,
                                         str_project_dir, f_is_hg_38):
        """

        * args_call: Arguments for the pipeline
                   : Dict
        * str_variants_file: Path to file to be annotated and filtered
                           : String path
        * f_is_hg_38: Indicates if the reference is Hg38, Hg19, or something else (in which CRAVAT will not run)
                    : True (HG38), False (HG19), None (CRAVAT will not run)
        * return: List of commands
        """

        # TODO
        # ClinVAR, CADD for annotation?
        # TODO, Do all COMIC IDs represent pathogenic mutations, if not remove non-pathogenic variants with COSMIC IDS

        # Commands for cancer filtering
        lcmd_cancer_filter = []

        # File to filter (may be annotated with cosmic or not so the name changes
        str_vcf_to_filter = str_variants_file

        # Files created
        str_vcf_base = os.path.join(str_project_dir, STR_MISC_DIR, os.path.basename(os.path.splitext(str_variants_file)[0]))
        str_cancer_mutations_unfiltered = os.path.join(str_project_dir, STR_MISC_DIR, C_STR_CANCER_ANNOTATED_VCF)
        str_cancer_mutations_filtered = self.func_switch_ext(str_vcf_base, "_cosmic_filtered.vcf")
        str_cravat_annotated_coding_vcf = self.func_switch_ext(str_vcf_base, "_cosmic_filtered_cravate_annotated_coding.vcf.gz")
        str_cravat_annotated_all_vcf = str_project_dir + os.path.sep + "annotated_min_filtered.vcf.gz"
        str_cravat_filtered_groom_vcf = str_project_dir + os.path.sep + C_STR_CANCER_VCF
        str_cancer_tab = str_project_dir + os.path.sep + C_STR_CANCER_TAB
        str_cravat_result_dir = self.func_switch_ext(str_vcf_base, "_cosmic_filtered_cravat_annotations.gz")
        str_extracted_cravat_dir = str_vcf_base + "_cosmic_filtered_cravat_annotations"
        str_cravat_detail_coding = os.path.join(str_extracted_cravat_dir, "Variant.Result.tsv")
        str_cravat_detail_noncoding = os.path.join(str_extracted_cravat_dir, "Variant_Non-coding.Result.tsv")
        str_cravat_detail_coding_updated = os.path.join(str_project_dir, STR_MISC_DIR, "Variant_result_updated.tsv")
        str_cravat_detail_noncoding_updated = os.path.join(str_project_dir, STR_MISC_DIR, "Variant_non_coding_result_updated.tsv")
        str_return_vcf = str_cancer_mutations_filtered

        # Index and bgzip vcf
        dict_csi = self.func_csi(str_vcf_to_filter)
        lcmd_cancer_filter.extend(dict_csi[INDEX_CMD])
        str_vcf_to_filter = dict_csi[INDEX_FILE]

        # Pull out and annotate Coding Cancer Mutations
        # Adding the following annotations from COSMIC
        # If the VCF does not have an annotation in COSMIC then it is dropped
        # GENE, COSMIC_ID, TISSUE, TUMOR, FATHMM, SOMATIC
        if args_call.str_cosmic_coding_vcf:
            # Annotate cancer variants with COSMIC
            str_cancer_annotation_command = " ".join(["bcftools", "annotate", "--output-type", "z",
                                                      "--annotations", args_call.str_cosmic_coding_vcf,
                                                      "--columns", "INFO/COSMIC_ID,INFO/TISSUE,INFO/TUMOR,INFO/FATHMM,INFO/SOMATIC",
                                                      "--output", str_cancer_mutations_unfiltered, str_vcf_to_filter])
            cmd_cosmic = Command.Command(str_cur_command = str_cancer_annotation_command,
                                          lstr_cur_dependencies = [args_call.str_cosmic_coding_vcf, str_vcf_to_filter],
                                          lstr_cur_products = [str_cancer_mutations_unfiltered])
            lcmd_cancer_filter.append(cmd_cosmic)
            str_vcf_to_filter = str_cancer_mutations_unfiltered

        # Filter variant before CRAVAT
        str_cancer_filter_command = " ".join(["filter_vcf_for_cancer.py",
                                              str_vcf_to_filter,
                                              str_cancer_mutations_filtered])
        cmd_cancer_filter = Command.Command(str_cur_command = str_cancer_filter_command,
                                             lstr_cur_dependencies = [str_vcf_to_filter],
                                             lstr_cur_products = [str_cancer_mutations_filtered])
        lcmd_cancer_filter.append(cmd_cancer_filter)

        # Annotate non-common with CRAVAT
        str_cmd_make_cravat_tab = " ".join(["java -jar gatk-package-4.0.1.2-local.jar", "VariantsToTable", "-R", args_call.str_genome_fa,"-V", str_cravat_filtered_groom_vcf,
                                           "-F", "CHROM", "-F", "POS", "-F", "REF", "-F", "ALT", "-F", "GENE",
                                           "-F", "DP", "-F", "QUAL", "-F", "MQ",
                                           "-F", "SAO", "-F", "NSF", "-F", "NSM", "-F", "NSN", "-F", "TUMOR", "-F", "TISSUE",
                                           "-F", "COSMIC_ID", "-F", "KGPROD", "-F", "RS", "-F", "PMC"])
        str_pred_filtered_vcf=str_cancer_mutations_filtered
        if (not f_is_hg_38 is None) and (not args_call.f_skip_cravat):
            # Update the output target vcf file given these steps are ran.
            str_return_vcf = str_cravat_filtered_groom_vcf
            str_pred_filtered_vcf=self.func_switch_ext(str_vcf_base, "_cosmic_filtered_cravate_annotated_filtered.vcf")
            str_cmd_make_cravat_tab = " ".join([str_cmd_make_cravat_tab,
                                                "-F", "CHASM_PVALUE",
                                                "-F", "CHASM_FDR",
                                                "-F", "VEST_PVALUE",
                                                "-F", "VEST_FDR"])

            str_cravat_result_dir_zip = str_cravat_result_dir + ".zip"
            lstr_hg_38 = ["--is_hg19"] if not f_is_hg_38 else []
            str_cravat_cmd = " ".join(["annotate_with_cravat.py", "--classifier", args_call.str_cravat_classifier] + lstr_hg_38 +
                                      ["--email", args_call.str_email_contact, "--max_attempts", str(I_CRAVAT_ATTEMPTS),
                                      "--wait", str(I_CRAVAT_WAIT), str_cancer_mutations_filtered, str_cravat_result_dir])
            cmd_cravat = Command.Command(str_cur_command = str_cravat_cmd,
                                      lstr_cur_dependencies = [str_cancer_mutations_filtered],
                                      lstr_cur_products = [str_cravat_result_dir_zip])
            lcmd_cancer_filter.append(cmd_cravat)

            ## Unzip
            str_unzip_cravat_cmd = " ".join(["unzip", "-d", str_extracted_cravat_dir, str_cravat_result_dir_zip])
            cmd_unzip_cravat = Command.Command(str_cur_command = str_unzip_cravat_cmd,
                                               lstr_cur_dependencies = [str_cravat_result_dir_zip],
                                               lstr_cur_products = [str_extracted_cravat_dir])
            lcmd_cancer_filter.append(cmd_unzip_cravat)

            # MV files needed from the CRAVAT dir to the current working dir.
            str_coding_variant_result = str_extracted_cravat_dir+os.path.sep+"*"+os.path.sep+"Variant.Result.tsv"
            str_non_coding_variant_result = str_extracted_cravat_dir+os.path.sep+"*"+os.path.sep+"Variant_Non-coding.Result.tsv"
            str_move_cravate_files = " ".join(["bash -c \"cp", "{"+str_coding_variant_result+","+str_non_coding_variant_result+"}", str_extracted_cravat_dir+"\""])
            cmd_mv_cravat = Command.Command(str_cur_command = str_move_cravate_files,
                                           lstr_cur_dependencies = [str_extracted_cravat_dir],
                                           lstr_cur_products = [str_cravat_detail_noncoding, str_cravat_detail_coding])
            lcmd_cancer_filter.append(cmd_mv_cravat)

            # Groom CRAVAT output tab for it does not upset BCFtools.
            str_groom_cravat_tab_coding = " ".join(["groom_cravat_annotation.py", str_cravat_detail_coding, str_cravat_detail_coding_updated])
            str_groom_cravat_tab_non_coding = " ".join(["groom_cravat_annotation.py", str_cravat_detail_noncoding, str_cravat_detail_noncoding_updated])
            cmd_groom_cravat_tab_coding = Command.Command(str_cur_command=str_groom_cravat_tab_coding,
                                                         lstr_cur_dependencies=[str_extracted_cravat_dir, str_cravat_detail_coding],
                                                         lstr_cur_products=[str_cravat_detail_coding_updated])
            cmd_groom_cravat_tab_noncoding = Command.Command(str_cur_command=str_groom_cravat_tab_non_coding,
                                                         lstr_cur_dependencies=[str_extracted_cravat_dir, str_cravat_detail_noncoding],
                                                         lstr_cur_products=[str_cravat_detail_noncoding_updated])
            lcmd_cancer_filter.extend([cmd_groom_cravat_tab_coding, cmd_groom_cravat_tab_noncoding])

            # Tabix index the CRAVAT tsv files
            dict_tabix = self.func_tabix(str_cravat_detail_coding_updated, str_output_dir = os.path.join(str_project_dir, STR_MISC_DIR), str_tabix = "-s 1 -b 2 -e 2 -S 12")
            lcmd_cancer_filter.extend(dict_tabix[INDEX_CMD])
            str_cravat_detail_coding_updated = str_cravat_detail_coding_updated +".gz"
            dict_tabix = self.func_tabix(str_cravat_detail_noncoding_updated, str_output_dir = os.path.join(str_project_dir, STR_MISC_DIR), str_tabix = "-s 1 -b 2 -e 2 -S 12")
            lcmd_cancer_filter.extend(dict_tabix[INDEX_CMD])
            str_cravat_detail_noncoding_updated = str_cravat_detail_noncoding_updated +".gz"

            ## Annotate and VCF file with TAB data.
            ## CRAVAT gives both Coding and none coding Variants results.
            ## For now, including both and not excluding noncoding.
            str_annotate_coding = " ".join(["bcftools",
                                            "annotate",
                                            "--annotations",
                                            str_cravat_detail_coding_updated,
                                            "-h", args_call.str_cravat_headers,
                                            "--columns",
                                            "\"CHROM,POS,CHASM_PVALUE,CHASM_FDR,VEST_PVALUE,VEST_FDR\"",
                                            "--output-type",
                                            "z",
                                            "--output",
                                            str_cravat_annotated_coding_vcf,
                                            str_cancer_mutations_filtered])
            cmd_annotate_with_cravat_coding = Command.Command(str_cur_command = str_annotate_coding,
                                                      lstr_cur_dependencies = [str_cravat_detail_coding_updated, str_cancer_mutations_filtered],
                                                      lstr_cur_products = [ str_cravat_annotated_coding_vcf])
            lcmd_cancer_filter.append(cmd_annotate_with_cravat_coding)
            str_annotate_noncoding = " ".join(["bcftools",
                                               "annotate",
                                               "--annotations",
                                               str_cravat_detail_noncoding_updated,
                                               "-h",
                                               args_call.str_cravat_headers,
                                               "--columns",
                                               "\"CHROM,POS,CHASM_PVALUE,CHASM_FDR,VEST_PVALUE,VEST_FDR\"",
                                               "--output-type",
                                               "z",
                                               "--output",
                                               str_cravat_annotated_all_vcf,
                                               str_cravat_annotated_coding_vcf])
            cmd_annotate_with_cravat_noncoding = Command.Command(str_cur_command = str_annotate_noncoding,
                                                      lstr_cur_dependencies = [str_cravat_detail_noncoding_updated, str_cravat_annotated_coding_vcf],
                                                      lstr_cur_products = [str_cravat_annotated_all_vcf])
            lcmd_cancer_filter.append(cmd_annotate_with_cravat_noncoding)

            # Filter based on Predictions
            str_cmd_filter_predictions = " ".join(["filter_vcf_for_predictions.py",
                                                   str_cravat_annotated_all_vcf,
                                                   str_pred_filtered_vcf])
            cmd_filter_predictions = Command.Command(str_cur_command = str_cmd_filter_predictions,
                                                     lstr_cur_dependencies = [str_cravat_annotated_all_vcf],
                                                     lstr_cur_products = [str_pred_filtered_vcf])
            lcmd_cancer_filter.append(cmd_filter_predictions)
            cmd_filter_predictions.func_set_dependency_clean_level([str_cravat_annotated_all_vcf], Command.CLEAN_NEVER)

        # Groom before filter
        str_cmd_groom_cancer_filtered = " ".join(["groom_vcf.py",
                                                  str_pred_filtered_vcf,
                                                  str_cravat_filtered_groom_vcf])
        cmd_groom_cancer_filtered = Command.Command(str_cur_command = str_cmd_groom_cancer_filtered,
                                                    lstr_cur_dependencies = [str_pred_filtered_vcf],
                                                    lstr_cur_products = [str_cravat_filtered_groom_vcf])
        lcmd_cancer_filter.append(cmd_groom_cancer_filtered)

        # Convert filtered VCF file to tab file.
        str_cmd_make_cravat_tab = " ".join([str_cmd_make_cravat_tab,
                                            "--lenient",
                                            "-O", str_cancer_tab])
        cmd_cravat_table = Command.Command(str_cur_command = str_cmd_make_cravat_tab,
                                           lstr_cur_dependencies = [str_cravat_filtered_groom_vcf],
                                           lstr_cur_products = [str_cancer_tab])
        cmd_cravat_table.func_set_dependency_clean_level([str_cravat_filtered_groom_vcf], Command.CLEAN_NEVER)
        lcmd_cancer_filter.append(cmd_cravat_table)
        return {INDEX_CMD:lcmd_cancer_filter,
                INDEX_FILE:str_cancer_tab}

    def func_make_commands(self, args_parsed, cur_pipeline):
        """
        Runs the pipeline. This is placed in a function so that multiple
        scripts with different arguments requirements can run it.
        For instance, one script managing the complete pipeline requiring
        many more arguments while having a simple script running only the
        first indexing step to make a global index useable in all subsequent
        runs associated with the reference genome generating the index.
        * args_parsed : Arguments
                    : Arguments used to run pipeline.

        """
        # Reset args_parsed if validating with DNA-seq data
        ### Turning off DNA Seq functionality
        args_parsed.f_validate_by_dnaseq = False
        #if args_parsed.f_validate_by_dnaseq:
        #    args_parsed.str_alignment_mode = STR_DNASEQ_VALIDATION
        #    args_parsed.str_variant_call_mode = STR_DNASEQ_VALIDATION
        #    args_parsed.f_recalibrate_sam = True

        # Fastq files or bam files must be given or index only should be true
        if(not (args_parsed.str_sample_file_left_fq
                and args_parsed.str_sample_file_right_fq)
           and not args_parsed.str_bam_file):
           str_error = "".join(["RNASEQ MUTATION PIPELINE, please make",
                                "sure to inclue a bam file or paired",
                                "fastq files unless running in index",
                                "only mode (which does no mutation calling)."])
           cur_pipeline.logr_logger.error(str_error)
           exit(7)

        # Constants
        # If not using a premade index and indexing only, do not update the
        # name of the index dir with the sample name. Does not need to be
        # unique. Otherwise update the name of the index directory so it is
        # unique, in case the pipeline is ran for multiple samples at once.
        str_sample_postfix = os.path.splitext(os.path.basename(args_parsed.str_sample_file_left_fq if args_parsed.str_sample_file_left_fq else args_parsed.str_bam_file))[0]
        str_sample_postfix = str_sample_postfix.replace(".","_")

        # Make the gtf files required for GSNAP flows
        #if args_parsed.str_alignment_mode in [STR_ALIGN_GSNP]:
        #    if not args_parsed.str_gtf_file_path:
        #        print("".join(["GTF file is required when using the ",
        #                       args_parsed.str_alignment_mode,
        #                       " alignment method. Please provide and try again."]))
        #        exit(5)

        # Vary alignment depending on arguments
        dict_align_funcs = {STR_ALIGN_STAR_LIMITED:self.func_do_star_alignment,
                            STR_ALIGN_STAR:self.func_do_star_alignment,
                            STR_DNASEQ_VALIDATION:self.func_do_BWA_alignment}

        # Vary variant calling depending on arguments
        dict_calling_funcs = {STR_DNASEQ_VALIDATION:self.func_call_dnaseq_like_rnaseq,
                              STR_VARIANT_GATK:self.func_do_variant_calling_gatk,
                              STR_VARIANT_SAMTOOLS:self.func_do_variant_calling_samtools,
                              STR_VARIANT_NONE:self.func_do_variant_calling_none }

        # Vary variant filtration
        dict_filtering_funcs = {STR_FILTERING_BCFTOOLS:self.func_do_filtering_bcftools,
                                STR_FILTERING_GATK:self.func_do_filtering_gatk,
                                STR_FILTERING_NONE:self.func_do_filtering_none}

        # If the output directory is not given,
        # get the file base from a sample file
        if args_parsed.f_wdl_run:
            args_parsed.str_out_dir = ""
        else:
            if not args_parsed.str_out_dir:
                args_parsed.str_out_dir = str_sample_postfix

        # Make sure the output directory is absolute
        args_parsed.str_out_dir = os.path.abspath(args_parsed.str_out_dir)

        # Base outputs on the sample file unless an output directory is given
        # Directories
        str_misc_dir = os.path.join(args_parsed.str_out_dir, STR_MISC_DIR)

        # Make pipeline object and indicate Log file
        #pline_cur = Pipeline.Pipeline(str_name = "rnaseq_mutation",
        #                              str_log_to_file = args_parsed.str_log_file,
        #                              str_update_source_path = args_parsed.str_update_classpath if hasattr(args_parsed, "str_update_classpath") else None)

        # Start commands
        lcmd_commands = []
        dict_align_info = {}

        if not args_parsed.str_bam_file:
            # If a bam file is given, ignore alignment and use the bam.
            # Handle indexing and alignment
            # Vary handling based on alignment type
            dict_align_info = dict_align_funcs[args_parsed.str_alignment_mode](args_call=args_parsed,
                                                                               str_unique_id=str_sample_postfix)

            # Run commands but indexing only
            lcmd_commands.extend(dict_align_info[INDEX_CMD])

        # Alignment method is previously used for indexing but at this point,
        # if a bam is given, the pipeline ignores alignment method and uses the bam
        if args_parsed.str_bam_file:
            dict_align_info = {INDEX_FILE:args_parsed.str_bam_file,
                               INDEX_FOLDER:args_parsed.str_bam_file}

        # Make directories
        #if not pline_cur.func_mkdirs([str_misc_dir]):
        #    exit(3)

        # If making depth files
        if args_parsed.f_calculate_base_coverage and (args_parsed.str_variant_call_mode == STR_VARIANT_NONE):
            # Create depth file
    #        str_depth_compressed_file = os.path.basename(self.func_switch_ext(dict_align_info[INDEX_FILE], "_depth.gz"))
            str_depth_compressed_file = os.path.basename(self.func_switch_ext(dict_align_info[INDEX_FILE], "depth"))
            str_depth_compressed_file = os.path.join(args_parsed.str_out_dir, str_depth_compressed_file)
    #        lcmd_commands.append(Command.Command(str_cur_command = "samtools depth " + dict_align_info[INDEX_FILE] + " | gzip > " + str_depth_compressed_file,
            lcmd_commands.append(Command.Command(str_cur_command = "samtools depth " + dict_align_info[INDEX_FILE] + " > " + str_depth_compressed_file,
                                                   lstr_cur_dependencies = [dict_align_info[INDEX_FILE]],
                                                   lstr_cur_products = [str_depth_compressed_file]))

        # If variant calling is occuring
        if(args_parsed.str_variant_call_mode.lower() != STR_VARIANT_NONE.lower()):

            # Currently edited VCF file
            str_annotated_vcf_file = ""

            str_json_inspector_file = args_parsed.str_out_dir + os.path.sep + C_STR_MUTATION_INSPECTOR

            # Add variant calling commands
            func_calling = dict_calling_funcs[args_parsed.str_variant_call_mode]
            dict_ret_variant_calling = func_calling(args_call=args_parsed,
                                                    str_align_file=dict_align_info[INDEX_FILE],
                                                    str_unique_id=str_sample_postfix,
                                                    str_project_dir=args_parsed.str_out_dir,
                                                    str_tmp_dir=str_misc_dir,
                                                    lstr_dependencies=[dict_align_info[INDEX_FILE],
                                                                       dict_align_info[INDEX_FILE]+".bai"],
                                                    logr_cur=cur_pipeline.logr_logger)
            lcmd_commands.extend(dict_ret_variant_calling[INDEX_CMD])
            str_bam_called_from = dict_ret_variant_calling[INDEX_BAM]
            str_bai_called_from = dict_ret_variant_calling[INDEX_INDEX]
            str_annotated_vcf_file = dict_ret_variant_calling[INDEX_FILE]

            # Add variant filtering
            func_filtering = dict_filtering_funcs[args_parsed.str_variant_filter_mode]
            dict_ret_variant_filtration = func_filtering(args_call=args_parsed,
                                                         str_variants_file=dict_ret_variant_calling[INDEX_FILE],
                                                         lstr_dependencies=[dict_ret_variant_calling[INDEX_FILE]],
                                                         logr_cur=cur_pipeline.logr_logger)
            if len(dict_ret_variant_filtration[INDEX_CMD]):
                lcmd_commands.extend(dict_ret_variant_filtration[INDEX_CMD])
                str_annotated_vcf_file = dict_ret_variant_filtration[INDEX_FILE]

            # If using the DNASEQ mode then stop,the rest is for RNASEQ.
            if not args_parsed.f_validate_by_dnaseq:

                # Clean up VCF file after variant caller
                str_clean_vcf = os.path.join(str_misc_dir,
                                             self.func_switch_ext(os.path.basename(str_annotated_vcf_file),
                                             "_clean.vcf"))
                str_clean_vcf_cmd = " ".join(["groom_vcf.py",
                                              str_annotated_vcf_file,
                                              str_clean_vcf])
                cmd_clean_vcf = Command.Command(str_cur_command=str_clean_vcf_cmd,
                                             lstr_cur_dependencies=[str_annotated_vcf_file],
                                             lstr_cur_products=[str_clean_vcf])
                cmd_clean_vcf.func_set_dependency_clean_level([str_annotated_vcf_file],
                                                              Command.CLEAN_NEVER)
                str_annotated_vcf_file = str_clean_vcf
                lcmd_commands.append(cmd_clean_vcf)

    #            lcmd_commands.extend(self.func_plot_vcf(str_annotated_vcf_file)[INDEX_CMD])

                # Filter results to just SNPs
                str_snp_filtered_vcf = self.func_switch_ext(str_annotated_vcf_file, "_snp.vcf")
                str_cmd_filter_snps = " ".join(["reduce_vcf_to_snps.py", str_annotated_vcf_file, str_snp_filtered_vcf])
                cmd_snp_filter = Command.Command(str_cur_command = str_cmd_filter_snps,
                                              lstr_cur_dependencies = [str_annotated_vcf_file],
                                              lstr_cur_products = [str_snp_filtered_vcf])
                lcmd_commands.append(cmd_snp_filter)
                str_annotated_vcf_file = str_snp_filtered_vcf
    #            lcmd_commands.extend(self.func_plot_vcf(str_annotated_vcf_file)[INDEX_CMD])

                # Filter RNA Editing
                if args_parsed.str_darned_data or args_parsed.str_radar_data:
                    str_rna_edit_filtered_vcf = self.func_switch_ext(str_annotated_vcf_file, "_RNAedit.vcf")
                    lstr_cmd_rna_editing_filter = ["filter_snps_rna_editing.py"]
                    if args_parsed.str_darned_data:
                        lstr_cmd_rna_editing_filter.extend(["--darned", args_parsed.str_darned_data])
                    if args_parsed.str_radar_data:
                        lstr_cmd_rna_editing_filter.extend(["--radar", args_parsed.str_radar_data])
                    lstr_cmd_rna_editing_filter.extend([str_annotated_vcf_file, str_rna_edit_filtered_vcf])
                    str_cmd_rna_editing_filter = " ".join(lstr_cmd_rna_editing_filter)
                    cmd_rna_editing_filter = Command.Command(str_cur_command = str_cmd_rna_editing_filter,
                                                          lstr_cur_dependencies = [str_annotated_vcf_file],
                                                          lstr_cur_products = [str_rna_edit_filtered_vcf])
                    lcmd_commands.append(cmd_rna_editing_filter)
                    # Switch over the annotated VCF to this RNA-Edited annotated VCF
                    str_annotated_vcf_file = str_rna_edit_filtered_vcf
    #                lcmd_commands.extend(self.func_plot_vcf(str_annotated_vcf_file)[INDEX_CMD])

                # Tabix / gz file sample
                dict_sample_csi = self.func_csi(str_annotated_vcf_file,
                                           args_parsed.str_out_dir)
                str_annotated_vcf_file = dict_sample_csi[INDEX_FILE]
                str_csi_vcf_file = dict_sample_csi[INDEX_INDEX]
                lcmd_commands.extend(dict_sample_csi[INDEX_CMD])

                # Tabix / gz DBSNP
                dict_dbsnp_csi = self.func_csi(args_parsed.str_vcf_file,
                                          args_parsed.str_out_dir)
                lcmd_commands.extend(dict_dbsnp_csi[INDEX_CMD])
                str_compressed_dbsnp = args_parsed.str_vcf_file + ".gz"
                str_csi_dbsnp = dict_dbsnp_csi[INDEX_INDEX]

                # DBSNP annotation
                # Annotate combined sample vcf files
                # bcftools annotate --annotations str_dbsnp_vcf -c
                # PM variant is clinicall precious (clinical and pubmed cited)
                # NSF, NSM, NSN, COMMON, SAO, KGPROD, KGVALIDATED, MUT, WTD, VLD, RS, PMC
                str_dbsnp_annotated_vcf = self.func_switch_ext(str_annotated_vcf_file, "_dbsnp.vcf.gz")
                str_annotate_command = "".join(["bcftools", " annotate",
                                                 " --output-type", " z",
                                                 " --annotations ",
                                                 str_compressed_dbsnp,
                                                 " --columns ",
                                                 "INFO/COMMON,INFO/PM,INFO/NSF,",
                                                 "INFO/NSM,INFO/NSN,INFO/SAO,",
                                                 "INFO/KGPROD,INFO/KGValidated,",
                                                 "INFO/MUT,INFO/WTD,INFO/VLD,",
                                                 "INFO/RS,INFO/PMC", " --output ",
                                                 str_dbsnp_annotated_vcf,
                                                 " ", str_annotated_vcf_file])
                lcmd_commands.append(Command.Command(str_cur_command=str_annotate_command,
                                               lstr_cur_dependencies=[str_compressed_dbsnp,
                                                                      str_csi_dbsnp,
                                                                      str_annotated_vcf_file,
                                                                      str_csi_vcf_file],
                                               lstr_cur_products=[str_dbsnp_annotated_vcf]))
                str_annotated_vcf_file = str_dbsnp_annotated_vcf

                # SNPeff java -jar /seq/regev_genome_portal/SOFTWARE/snpEff/snpEff.jar -nostats -noLof -no-downstream -no-upstream hg19 variants.vcf > new.vcf
                str_snp_eff_annotated = self.func_switch_ext(str_annotated_vcf_file, "_snpeff.vcf")
                str_snp_eff_cmd = " ".join(["bgzip -cd", str_annotated_vcf_file,
                                            "|", "java -jar snpEff.jar -nostats",
                                            "-noLof -no-downstream -no-upstream",
                                            "hg19", ">", str_snp_eff_annotated])
                lcmd_commands.append(Command.Command(str_cur_command = str_snp_eff_cmd,
                                                   lstr_cur_dependencies = [str_annotated_vcf_file],
                                                   lstr_cur_products = [str_snp_eff_annotated]))
                str_annotated_vcf_file = str_snp_eff_annotated

                # Update the SNPeff style annotations to the simple info column feature style
                str_snp_eff_updated_file = self.func_switch_ext(str_annotated_vcf_file, "_updated.vcf")
                str_snp_eff_update_cmd = " ".join(["update_snpeff_annotations.py", str_annotated_vcf_file, str_snp_eff_updated_file])
                lcmd_commands.append(Command.Command(str_cur_command = str_snp_eff_update_cmd,
                                                   lstr_cur_dependencies = [str_annotated_vcf_file],
                                                   lstr_cur_products = [str_snp_eff_updated_file]))
                str_annotated_vcf_file = str_snp_eff_updated_file

                # Perform cancer filtering
                f_cravat_hg38 = None
                if args_parsed.str_email_contact is None or (not args_parsed.f_hg_19 and not args_parsed.f_hg_38):
                    cur_pipeline.logr_logger.warning("CRAVAT analysis will not be ran. Please make sure to provide an email and indicate if hg38 or hg19 is being used.")
                elif args_parsed.f_hg_38:
                    f_cravat_hg38 = True
                elif args_parsed.f_hg_19:
                    f_cravat_hg38 = False
                cmd_filter_cancer = self.func_do_variant_filtering_cancer(args_call=args_parsed,
                                                                          str_variants_file=str_annotated_vcf_file,
                                                                          str_project_dir=args_parsed.str_out_dir,
                                                                          f_is_hg_38=f_cravat_hg38)
                str_cancer_tab = cmd_filter_cancer[INDEX_FILE]
                lcmd_commands.extend(cmd_filter_cancer[INDEX_CMD])

                # Make JSON file for the inspector
                if args_parsed.str_bed:
                    str_cmd_json_inspector = " ".join(["make_mutation_inspector_json.py",
                                                "--sample", args_parsed.str_out_dir,
                                                "--tab", str_cancer_tab,
                                                "--bam", str_bam_called_from,
                                                "--bam_index", str_bai_called_from,
                                                "--bed", args_parsed.str_bed,
                                                "--bed_index", args_parsed.str_bed + ".idx",
                                                str_json_inspector_file])
                    cmd_json_inspector = Command.Command(str_cur_command = str_cmd_json_inspector,
                                                 lstr_cur_dependencies = [str_cancer_tab, str_bam_called_from, str_bai_called_from],
                                                 lstr_cur_products = [str_json_inspector_file])
                    cmd_json_inspector.func_set_dependency_clean_level([str_cancer_tab, str_bam_called_from, str_bai_called_from], Command.CLEAN_NEVER)
                    lcmd_commands.append(cmd_json_inspector)

                    # Copy bed to output to make it an output for Galaxy and allow it to be used in the inspector.
                    str_copied_bed = os.path.join(args_parsed.str_out_dir, os.path.basename(args_parsed.str_bed))
                    str_cmd_copy_bed = " ".join(["cp", args_parsed.str_bed, str_copied_bed])
                    lcmd_commands.append(Command.Command(str_cur_command = str_cmd_copy_bed,
                                                     lstr_cur_dependencies = [args_parsed.str_bed],
                                                     lstr_cur_products = [str_copied_bed]))

        # Run commands including variant calling
        return(lcmd_commands)

    def func_gz(self, str_vcf, str_output_dir = ""):
      """
          Creates a bcftools index (vcf index) for the given vcf file.
          If it is not gzipped, the directory it is in is checked for a gz file.
          If the gz file does not exist the file is gzipped.
      """
      lcmd_index = []

      # Check extension
      if not os.path.splitext(str_vcf)[1] == ".gz":

        # Check if the gz file exists
        if os.path.exists(str_vcf + ".gz"):
          str_vcf = str_vcf + ".gz"
        else:
          # GZ files
          str_gz = str_vcf + ".gz"
          if str_output_dir:
            str_gz = os.path.join(str_output_dir, os.path.basename(str_gz))
          str_cmd_gz = " ".join(["bgzip -c ", str_vcf, ">", str_gz])
          cmd_gz = Command.Command(str_cur_command = str_cmd_gz,
                                           lstr_cur_dependencies = [str_vcf],
                                           lstr_cur_products = [str_gz])
          str_vcf = str_gz
          lcmd_index.append(cmd_gz)
      return({ INDEX_CMD: lcmd_index, INDEX_FILE: str_vcf })


    def func_csi(self, str_vcf, str_output_dir = ""):
      """
          Creates a bcftools index (vcf index) for the given vcf file.
          If it is not gzipped, the directory it is in is checked for a gz file.
          If the gz file does not exist the file is gzipped.
          A csi file is also made if it does not exist using the gz file.
      """
      lcmd_index = []

      # Gzip
      dict_gz = self.func_gz(str_vcf, str_output_dir)
      str_vcf = dict_gz[INDEX_FILE]
      if dict_gz[INDEX_CMD]:
        lcmd_index.extend(dict_gz[INDEX_CMD])

      str_index_file = str_vcf + ".csi"
      if not os.path.exists(str_index_file):
        # Create index for the VCF file
        if str_output_dir:
          str_index_file = os.path.join(str_output_dir, os.path.basename(str_index_file))
        str_cmd_index_vcf = " ".join(["bcftools index", str_vcf])
        cmd_index_vcf = Command.Command(str_cur_command = str_cmd_index_vcf,
                                         lstr_cur_dependencies = [str_vcf],
                                         lstr_cur_products = [str_index_file])
        lcmd_index.append(cmd_index_vcf)
      return({INDEX_CMD:lcmd_index,
              INDEX_FILE:str_vcf,
              INDEX_INDEX:str_index_file})

    def func_tabix(self, str_vcf, str_output_dir = "", str_tabix = ""):
      """
          Creates a tbi (vcf index) for the given vcf file.
          If it is not gzipped, the directroy is checked for the gz file.
          If the gz file does not exist, the file is gzipped.
          The tbi file is then made from the gz file.
      """
      lcmd_tabix = []

      # Gzip
      dict_gz = self.func_gz(str_vcf, str_output_dir)
      str_vcf = dict_gz[INDEX_FILE]
      if dict_gz[INDEX_CMD]:
        lcmd_tabix.extend(dict_gz[INDEX_CMD])

      if not os.path.exists(str_vcf + ".tbi"):
        # Create index for the VCF file
        str_tbi = str_vcf + ".tbi"
        if str_output_dir:
          str_tbi = os.path.join(str_output_dir, os.path.basename(str_tbi))
        str_cmd_tabix_vcf = " ".join(["tabix -f",str_tabix, str_vcf])
        cmd_tabix_vcf = Command.Command(str_cur_command = str_cmd_tabix_vcf,
                                         lstr_cur_dependencies = [str_vcf],
                                         lstr_cur_products = [str_tbi])
        lcmd_tabix.append(cmd_tabix_vcf)
      return({ INDEX_CMD: lcmd_tabix, INDEX_FILE: str_vcf })

    def func_plot_vcf(self, str_vcf):
      lcmd_plot = []
      str_plot_location = os.path.join(os.path.dirname(str_vchk_stats), os.path.basename(str_vchk_stats) + "_plot")
      str_vchk_stats_command = " ".join(["bcftools", "stats", str_vcf, ">", str_vchk_stats])
      str_vchk_plot_command = " ".join(["plot-vcfstats", str_vchk_stats, "-p", str_plot_location + os.path.sep])
      lcmd_plot.append(Command.Command(str_cur_command = str_vchk_stats_command,
                                          lstr_cur_dependencies = [str_vcf],
                                          lstr_cur_products = [str_vchk_stats]))
      lcmd_plot.append(Command.Command(str_cur_command = str_vchk_plot_command,
                                          lstr_cur_dependencies = [str_vchk_stats],
                                          lstr_cur_products = [str_plot_location]))
      return({ INDEX_CMD: lcmd_plot, INDEX_FILE: str_plot_location })

    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """
        # Parse arguments
        arg_raw.prog = "rnaseq_mutation_pipeline.py"
        arg_raw.description = "Variant calling using RNASeq NGS sequencing"
        arg_raw.add_argument("-l", "--left", metavar = "Left_sample_file", dest = "str_sample_file_left_fq", required = False, help = "Path to one of the two paired RNAseq samples (left)")
        arg_raw.add_argument("-r", "--right", metavar = "Right_sample_file", dest = "str_sample_file_right_fq", required = False, help = "Path to one of the two paired RNAseq samples (right)")
        arg_raw.add_argument("-n", "--threads", metavar = "Process_threads", dest = "i_number_threads", type = int, default = 1, help = "The number of threads to use for multi-threaded steps.")
        arg_raw.add_argument("--wdl_compatible_run", dest = "f_wdl_run", action="store_true", help = "Cromwell/WDL requires execution to happen relative to an output directory of dynaically created without giving the directory path to the underlying tools/pipelines. This requires a pipeline to use relative paths which can be dangerous outside of Cromwell/WDL. This will ignore any output directory specified in the command line and force the output to be relative paths. DO NOT USE outside of Cromwell/WDL.")

        # Run modes
        args_group_run = arg_raw.add_argument_group("Run Mode", "Associated in running different modes of the pipeline.")
        args_group_run.add_argument("--bam", metavar = "bam_file", dest = "str_bam_file", default = None, help = "Sample file in the form of a bam, if this is given NO alignment will be performed; the alignment mode command line will be ignored; let and right sample files will be ignored. Normal pipeline processing will pick up directly after alignment in the pipeline with the supplied bam.")
        args_group_run.add_argument("-d", "--alignment_mode", metavar = "Alignment_mode", dest = "str_alignment_mode", default = STR_ALIGN_STAR, choices = LSTR_ALIGN_CHOICES, help = "Specifies the alignment and indexing algorithm to use.")
        args_group_run.add_argument("-e", "--variant_call_mode", metavar = "Call_mode", dest = "str_variant_call_mode", default = STR_VARIANT_GATK, choices = LSTR_VARIANT_CALLING_CHOICES, help = "Specifies the variant calling method to use.")
        args_group_run.add_argument("--variant_filtering_mode", metavar = "Filter_mode", dest = "str_variant_filter_mode", default = STR_FILTERING_DEFAULT, choices = LSTR_VARIANT_FILTERING_CHOICES, help = "Specifies the variant filtering method.")
        args_group_run.add_argument("--base_depth", dest = "f_calculate_base_coverage", default = False, action = "store_true", help = "Calculates the base coverage per base.")
        args_group_run.add_argument("-i", "--index", metavar = "Use_premade_index", dest = "str_initial_index", default = None, help = "The initial index is made only from the reference genome and can be shared. If premade, supply a path here to the index directory so that it is not rebuilt for every alignment. Please provide the full path.")
        #args_group_run.add_argument("--validate_dnaseq", dest = "f_validate_by_dnaseq", default = False, action = "store_true", help = "Used for development only. Should not be used with biological samples.")
        args_group_run.add_argument("-y", "--star_memory", metavar = "Star_memory", dest = "str_star_memory_limit", default = None, help = "Memory limit for star index. This should be used to increase memory if needed. Reducing memory consumption should be performed with the STAR Limited mod.")
        # GATK associated
        args_group_gatk = arg_raw.add_argument_group("GATK", "Associated with or controlling GATK tools.")
        args_group_gatk.add_argument("-a", "--realign", dest = "f_stop_optional_realignment", default = False, action = "store_true", help = "Turns off optional indel realignment step.")
        args_group_gatk.add_argument("-j", "--recalibrate_sam", dest = "f_recalibrate_sam", default = True, action="store_false", help = "If used, turns off gatk recalibration of bam files before samtools variant calling.")
        args_group_gatk.add_argument("-p", "--plot", dest = "f_optional_recalibration_plot", default = True, action = "store_false", help = "Turns off plotting recalibration of alignments.")
        arg_raw.add_argument("-s", "--sequencing_platform", metavar = "Sequencing Platform", dest = "str_sequencing_platform", default = "ILLUMINA", choices = LSTR_SEQ_CHOICES, help = "The sequencing platform used to generate the samples choices include " + " ".join(LSTR_SEQ_CHOICES) + ".")
        # Resources
        args_group_resources = arg_raw.add_argument_group("Resources", "Associated with resources for the pipelines.")
        args_group_resources.add_argument("--cosmic_vcf", metavar="cosmic_reference_vcf", dest="str_cosmic_coding_vcf", default=None, action="store", help="Coding Cosmic Mutation VCF annotated with Phenotype Information.")
        args_group_resources.add_argument("-f", "--reference", metavar = "Reference_genome", dest = "str_genome_fa", required = True, help = "Path to the reference genome to use in the analysis pipeline.")
        args_group_resources.add_argument("-k", "--gtf", metavar = "Reference GTF", dest = "str_gtf_file_path", default = None, help = "GTF file for reference genome.")
        args_group_resources.add_argument("-w", "--vcf", metavar = "Variant_calling_file_for_the_reference_genome", dest = "str_vcf_file", default = None, help = "Variant calling file for the reference genome.")
        args_group_resources.add_argument("--vcf_snps", metavar="Reference_VCF_of_SNPs", dest="str_snp_vcf", default=None, help="The reference VCF including only SNP entries, if not given, will be generated from the file given with --vcf")
        args_group_resources.add_argument("--darned", metavar = "Darned_data", dest = "str_darned_data", default = None, help = "Darned data for RNA editing removal, if included will be used for RNA editing removal.")
        args_group_resources.add_argument("--radar", metavar = "Radar_data", dest = "str_radar_data", default = None, help = "Radar data for RNA editing removal, if included will be used for RNA editing removal.")
        args_group_resources.add_argument("--bed", metavar = "Reference BED file", dest = "str_bed", default = None, help = "Bed file for reference genome, required only if making the mutation inspector json. If given the json file will be made. Please make sure the bed file is indexed and that bed.idx file is in the same folder with the same file basename as the related bed file.")
        # Cravat associated
        args_group_cravat = arg_raw.add_argument_group("CRAVAT", "Associated with CRAVAT prioritization of variant calls.")
        args_group_cravat.add_argument("--cravat_annotation_header", metavar = "cravat_headers", dest = "str_cravat_headers", default = None, help = "Headers for each CRAVAT feature annotated to the VCF file (used in BCFtools).")
        args_group_cravat.add_argument("--skip_cravat", dest = "f_skip_cravat", action="store_true", default = False, help = "Skips CRAVAT services.")
        args_group_cravat.add_argument("--tissue_type", metavar = "cravat_tissue", dest = "str_cravat_classifier", default = STR_CRAVAT_CLASSIFIER_DEFAULT, help = "Tissue type (used in CRAVAT variant prioritation). Supported classifiers can be found at http://www.cravat.us/help.jsp)")
        args_group_cravat.add_argument("--email", metavar = "email_contact", dest = "str_email_contact", default = None, help = "Email used to notify of errors associated with cravat.")
        group_hg = args_group_cravat.add_mutually_exclusive_group()
        group_hg.add_argument("--is_hg19", dest = "f_hg_19", action="store_true", help = "Indicates that Hg19 is being used.")
        group_hg.add_argument("--is_hg38", dest = "f_hg_38", action="store_true", help = "Indicates that Hg38 is being used.")
        return(arg_raw)

if __name__ == "__main__":

    # Need to rn, call the script
    RnaseqSnp().func_run_pipeline()
