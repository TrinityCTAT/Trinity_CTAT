#!/usr/bin/env python


__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2014"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"


import datetime
import os
import sciedpiper.Command as Command
import sciedpiper.Pipeline as Pipeline

# Constants
# Keys for alignment return ( dicts )
INDEX_CMD = "cmd"
INDEX_FILE = "out_file"
INDEX_FOLDER = "align_folder"

# Choices for platform
STR_ILLUMINA = "ILLUMINA"
LSTR_SEQ_CHOICES = [ STR_ILLUMINA, "SLX,SOLEXA", "SOLID,454", "COMPLETE", "PACBIO", "IONTORRENT", "CAPILLARY", "HELICOS" ]

# Choices for alignment
STR_ALIGN_GSNP = "GSNAP"
STR_ALIGN_STAR = "STAR"
STR_ALIGN_STAR_LIMITED = "LIMITED"
STR_ALIGN_TOPHAT = "TOPHAT"
LSTR_ALIGN_CHOICES = [ STR_ALIGN_STAR, STR_ALIGN_STAR_LIMITED ]

# Choices for variant calling
STR_VARIANT_GATK = "GATK"
STR_VARIANT_SAMTOOLS = "SAM"
LSTR_VARIANT_CALLING_CHOICES = [ STR_VARIANT_GATK ]

# This mode is used in validating the method in teh context fo DNA-seq
# It is not intended to be ran on biological samples for studies.
STR_DNASEQ_VALIDATION = "DNASEQ"


def func_do_star_alignment( args_call, str_unique_id, pline_cur, f_index_only = False ):
    """
    Manages the calls for star alignment.
    
    args_call : Arguments
              : Arguments used to run pipeline
              
    Return : List of commands
           : List of commands to run for star alignment
    """
    
    STR_ALIGN_1 = "_".join( [ "star_align_1", str_unique_id ] )
    STR_ALIGN_2 = "_".join( [ "star_align_2", str_unique_id ] )
    STR_INDEX_1 = "_".join( [ "star_index_1", str_unique_id ] )
    STR_INDEX_2 = "_".join( [ "star_index_2", str_unique_id ] )
    STR_STAR_GENOME_GENERATE = "genomeGenerate"
    STR_STAR_SPLICE_JUNCTION_FILE = "SJ.out.tab"
    
    # Update the limitGenomeGenerateRAM if more memory is requested
    lstr_index_memory_size = []
    if hasattr( args_call, "str_star_memory_limit" ):
        lstr_index_memory_size = [ "--limitGenomeGenerateRAM ", args_call.str_star_memory_limit ] if ( ( not ( args_call.str_star_memory_limit is None ) ) and 
                                                                                                   ( args_call.str_alignment_mode == STR_ALIGN_STAR ) ) else []
    
    lstr_limited_index_mode = [ "--limitGenomeGenerateRAM 15000000000 --genomeSAsparseD 2 --limitIObufferSize = 150000000" 
                                   ] if args_call.str_alignment_mode == STR_ALIGN_STAR_LIMITED else []
    lstr_limited_alignment_mode = [ "--genomeSAsparseD 2" ] if args_call.str_alignment_mode == STR_ALIGN_STAR_LIMITED else []

    str_align_dir_1 = os.path.join( args_call.str_file_base, STR_ALIGN_1 )
    str_align_dir_2 = os.path.join( args_call.str_file_base, STR_ALIGN_2 )
    str_index_dir_1 = os.path.join( args_call.str_file_base, STR_INDEX_1 )
    str_index_dir_1_change = os.path.join( "..", STR_INDEX_1 )
    str_index_dir_2 = os.path.join( args_call.str_file_base, STR_INDEX_2 )
    str_star_output_sam = os.path.join( str_align_dir_2, "Aligned.out.sam" )

    # Commands to build and return
    lcmd_commands = []

    # Make needed directories
    lstr_dirs_to_make = []
    # Handle premade index or build index scenarios
    if not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        lstr_dirs_to_make.append( str_index_dir_1 )
    else:
        STR_INDEX_1 = args_call.str_initial_index
        str_index_dir_1 = args_call.str_initial_index
        str_index_dir_1_change = args_call.str_initial_index
            
    # Handle not just indexing, also aligning
    if not f_index_only:
        lstr_dirs_to_make.extend( [ str_index_dir_2, str_align_dir_1, str_align_dir_2 ] )
    if not pline_cur.func_mkdirs( lstr_dirs_to_make ):
        exit( 3 )

    # If the premade index is not given then generate
    if not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "STAR --runMode", STR_STAR_GENOME_GENERATE ] +
                                                                         lstr_limited_index_mode + lstr_index_memory_size + [ "--genomeDir", str_index_dir_1, 
                                 "--genomeFastaFiles", args_call.str_genome_fa, "--runThreadN", str( args_call.i_number_threads ) ] ),
                                                lstr_cur_dependencies = [ args_call.str_genome_fa ],
                                                lstr_cur_products = [ str_index_dir_1 ] ) ] )

    # If aligning should occur as well
    if not f_index_only:

        # Update input files in case they are gzipped
        str_unzip_left, str_unzip_right = pline_cur.func_handle_gzip( [ os.path.join( "..", "..", args_call.str_sample_file_left_fq ), 
                                                                                        os.path.join( "..", "..", args_call.str_sample_file_right_fq ) ])
        # First star alignment
        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "cd", str_align_dir_1 ] ) ),
                            Command.Command( str_cur_command = " ".join( [ "STAR --genomeDir", str_index_dir_1_change ] +
                                                                         lstr_limited_alignment_mode + [ "--readFilesIn",
                                                                      str_unzip_left, str_unzip_right,
                                                                      "--runThreadN", str( args_call.i_number_threads ) ] ),
                                            lstr_cur_dependencies = [ args_call.str_sample_file_left_fq, 
                                                                 args_call.str_sample_file_right_fq, 
                                                                 str_index_dir_1 ],
                                            lstr_cur_products = [ str_align_dir_1 ] ) ] )

        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "cd", os.path.join( "..", ".." ) ] ) ),
                            # Second star index
                            Command.Command( str_cur_command = " ".join( [ "STAR --runMode", STR_STAR_GENOME_GENERATE ] +
                                                                         lstr_limited_index_mode + lstr_index_memory_size + [ "--genomeDir", str_index_dir_2,
                                                                      "--genomeFastaFiles", args_call.str_genome_fa, "--sjdbFileChrStartEnd", 
                                                                      os.path.join( str_align_dir_1, STR_STAR_SPLICE_JUNCTION_FILE ), 
                                                                      "--sjdbOverhang 75 --runThreadN", str( args_call.i_number_threads ) ] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, str_align_dir_1 ],
                                            lstr_cur_products = [ str_index_dir_2 ] ),
                            Command.Command( str_cur_command = " ".join( [ "cd", str_align_dir_2 ] ) ),
                            # Second star alignment
                            Command.Command( str_cur_command = " ".join( [ "STAR --genomeDir", os.path.join( "..", STR_INDEX_2 ) ] +
                                                                         lstr_limited_alignment_mode + [ "--readFilesIn",
                                                                      str_unzip_left, str_unzip_right,
                                                                      "--runThreadN", str( args_call.i_number_threads ) ] ),
                                            lstr_cur_dependencies = [ str_index_dir_2,
                                                                args_call.str_sample_file_left_fq,
                                                                args_call.str_sample_file_right_fq ],
                                            lstr_cur_products = [ str_align_dir_2 ] ),
                            Command.Command( str_cur_command = " ".join( [ "cd", os.path.join( "..", ".." ) ] ) ) ] )
    return { INDEX_CMD : lcmd_commands, INDEX_FILE : str_star_output_sam, INDEX_FOLDER : str_align_dir_2 }


def func_do_top_hat_alignment( args_call, str_unique_id, pline_cur, f_index_only = False ):
    """
    Manages the calls for the top hat alignment.
    
    args_call : Arguments
              : Arguments used to run pipeline
              
    Return : List of commands
           : List of commands to run for top hat alignment
    """
    
    # Holds the commands which will eventually be returned.
    lcmd_commands = []

    STR_INDEX = "_".join( [ "tophat_index", str_unique_id ] )
    STR_ALIGN = "_".join( [ "tophat_align", str_unique_id ] )
    STR_BAM_FILE = "accepted_hits.bam"
    str_index_dir = os.path.join( args_call.str_file_base, STR_INDEX )
    str_bt_index = os.path.join( str_index_dir, os.path.splitext( os.path.basename( args_call.str_genome_fa ) )[ 0 ] )
    str_align_dir = os.path.join( args_call.str_file_base, STR_ALIGN )
    str_output_file = os.path.join( str_align_dir, STR_BAM_FILE)

    # Make directories
    lstr_dirs_to_make = []
    # Handle premade index or build index scenarios
    if not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        lstr_dirs_to_make.append( str_index_dir )
    else:
        str_index_dir = args_call.str_initial_index
    # Handle not just indexing, also aligning
    if not f_index_only:
        lstr_dirs_to_make.append( str_align_dir )
    if not pline_cur.func_mkdirs( lstr_dirs_to_make ):
        exit( 5 )
    
    # Optionally make index
    if not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "bowtie2-build", args_call.str_genome_fa,
                                                                    str_bt_index ] ),
                                       lstr_cur_dependencies = [  args_call.str_genome_fa ],
                                       lstr_cur_products = [ str_index_dir ] ) ] )
        
    # Align
    if not f_index_only:

        # Update input files in case they are gzipped
        str_unzip_left, str_unzip_right = pline_cur.func_handle_gzip( [ os.path.join( "..", "..", args_call.str_sample_file_left_fq ), 
                                                                                        os.path.join( "..", "..", args_call.str_sample_file_right_fq ) ])

        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "tophat",
                                                                    "-o", str_align_dir,
                                                                    "-G", args_call.str_gtf_file_path,
                                                                    "--num-threads", str( args_call.i_number_threads ),
                                                                    os.path.join( str_index_dir, 
                                                                    os.path.splitext( os.path.basename( args_call.str_genome_fa ) )[0] ),
                                                                    str_unzip_left, str_unzip_right ] ),
                                       lstr_cur_dependencies = [ str_index_dir, args_call.str_gtf_file_path,
                                                                args_call.str_sample_file_left_fq,
                                                                args_call.str_sample_file_right_fq ],
                                       lstr_cur_products = [ str_align_dir ] ) ] )
    return { INDEX_CMD : lcmd_commands, INDEX_FILE : str_output_file , INDEX_FOLDER : str_align_dir }
    

def func_do_gsnp_alignment( args_call, str_unique_id, pline_cur, f_index_only = False ):
    """
    Manages the calls for gsnp alignment.
    
    args_call : Arguments
              : Arguments used to run pipeline
              
    Return : List of commands
           : List of commands to run for gsnp alignment
    """

    # Holds the commands which will eventually be returned.
    lcmd_commands = []

    STR_INDEX = "_".join( [ "gsnap_index", str_unique_id ] )
    STR_ALIGN = "_".join( [ "gsnap_align", str_unique_id ] )
    STR_BAM_FILE = "gsnap_align.sam"
    STR_IIT_STORE = "iit_store"
    STR_SPLICESITES = "gtf_splicesites"
    str_index_dir = os.path.join( args_call.str_file_base, STR_INDEX )
    str_genome_name = os.path.splitext( os.path.basename( args_call.str_genome_fa ) )[0]
    str_align_dir = os.path.join( args_call.str_file_base, STR_ALIGN )
    str_output_file = os.path.join( str_align_dir, STR_BAM_FILE)
    #TODO update, is sanger appropriate for all others?
    str_quality_protocol = "illumina" if args_call.str_sequencing_platform == STR_ILLUMINA else "sanger"

    # Make directories
    lstr_dirs_to_make = []
    # Handle premade index or build index scenarios
    if not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        lstr_dirs_to_make.append( str_index_dir )
    else:
        str_index_dir = args_call.str_initial_index
    # Handle not just indexing, also aligning
    if not f_index_only:
        lstr_dirs_to_make.append( str_align_dir )
    if not pline_cur.func_mkdirs( lstr_dirs_to_make ):
        exit( 5 )
    
    # Optionally make index
    if not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        # Directory made during the build process that has gene model mappings
        str_gmap_index = os.path.splitext( os.path.basename( str_genome_name ) )[ 0 ]
        str_gmap_map_dir = os.path.join( str_index_dir, str_gmap_index, ".".join( [ str_gmap_index, "maps" ] ) )
        str_splicesites_file = os.path.join( str_gmap_map_dir, ".".join( [ str_gmap_index, "splicesites" ] ) )
        str_splicesites_iit_file = ".".join( [ str_splicesites_file, "iit" ] )
        
        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "gmap_build", "-d", str_genome_name, 
                                                                    "-D", str_index_dir, args_call.str_genome_fa ] ),
                                       lstr_cur_dependencies = [  args_call.str_genome_fa ],
                                       lstr_cur_products = [ os.path.join( str_index_dir, str_genome_name ) ] ),
                                Command.Command( str_cur_command = " ".join( [ "cat", args_call.str_gtf_file_path, 
                                                                    "|", STR_SPLICESITES, ">", str_splicesites_file ] ),
                                                 lstr_cur_dependencies = [ args_call.str_gtf_file_path ],
                                                 lstr_cur_products = [ str_splicesites_file ] ),
                                 Command.Command( str_cur_command = " ".join( [ "cat", str_splicesites_file,
                                                                               "|", STR_IIT_STORE, "-o", str_splicesites_file ] ),
                                                  lstr_cur_dependencies = [ str_splicesites_file ],
                                                  lstr_cur_products = [ str_splicesites_iit_file ] ) ] )

    # Align
    if not f_index_only:

        # Update input files in case they are gzipped
        str_unzip_left, str_unzip_right = pline_cur.func_handle_gzip( [ os.path.join( "..", "..", args_call.str_sample_file_left_fq ), 
                                                                                        os.path.join( "..", "..", args_call.str_sample_file_right_fq ) ])   

        lcmd_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "cd", str_align_dir ] ) ),
                               Command.Command( str_cur_command = " ".join( [ "gsnap", "-d", str_genome_name,
                                                                    "-D", os.path.join( "..", "..", str_index_dir, str_genome_name ), 
                                                                    "-t", str( args_call.i_number_threads ),
                                                                    "-B 5 -a paired -N 1 -m 4 -M 1 -E 4 -n 100 --format=sam",
                                                                    "--gmap-mode=pairsearch,terminal,improve -O --quality-protocol="+str_quality_protocol,
                                                                    "-s", os.path.basename( str_splicesites_iit_file ),
                                                                    str_unzip_left, str_unzip_right, ">", "gsnap_align.sam" ] ),
                                       lstr_cur_dependencies = [ os.path.join( str_index_dir, str_genome_name ), 
                                                                args_call.str_sample_file_left_fq,
                                                                args_call.str_sample_file_right_fq ],
                                       lstr_cur_products = [ os.path.join( "..", "..", str_align_dir ) ] ),
                               Command.Command( str_cur_command = " ".join( [ "cd", os.path.join( "..", ".." ) ] ) ) ] )
    return { INDEX_CMD : lcmd_commands, INDEX_FILE : str_output_file , INDEX_FOLDER : str_align_dir }


def func_do_BWA_alignment( args_call, str_unique_id, pline_cur, f_index_only = False ):
    """
    Manages the calls for aligning DNA-seq data with BWA.
    Development use only, not intended to be used with studies.
    Used to validate technical properties of the RNA-seq pipeline.
    
    args_call : Arguments
              : Arguments used to run pipeline
              
    Return : List of commands
           : List of commands to run for BWA alignment
    """
    
    # Files
    str_left_file_key = os.path.basename( os.path.splitext( args.str_sample_file_left_fq )[ 0 ] )
    str_sam = os.path.join( args_call.str_file_base, ".".join( [ str_left_file_key, "sam" ] ) )
    str_bam = os.path.join( args_call.str_file_base, ".".join( [ str_left_file_key, "bam" ] ) )
    
    lcmd_dna_mapping_commands = []
    
    # Pre-processing
    ## Mapping and Dedupping
    ### BWA, make coordinate ordered bam
    #### Index reference
    if f_index_only or not hasattr( args_call, "str_initial_index" ) or not args_call.str_initial_index:
        cmd_bwa_index = Command.Command( str_cur_command = "".join( [ "bwa index -a bwtsw ", args.str_genome_fa ] ),
                                            lstr_cur_dependencies = [ args.str_genome_fa ],
                                            lstr_cur_products = [ args.str_genome_fa + ".amb",
                                                                 args.str_genome_fa + ".ann",
                                                                 args.str_genome_fa + ".bwt",
                                                                 args.str_genome_fa + ".pac",
                                                                 args.str_genome_fa + ".sa" ] )
        if f_index_only:
            return { INDEX_CMD: [ cmd_bwa_index ], INDEX_FOLDER: os.path.dirname( args.str_genome_fa ) }
        else:
            lcmd_dna_mapping_commands.append( cmd_bwa_index )

    # Align both samples
    # bwa sample ref.fasta fwd.sai rev.sai fwd.fq rev.fq > mydata.sam
    cmd_bwa_sam = Command.Command( str_cur_command = "".join( [ "bwa mem -M -R \"@RG\\tID:", str_unique_id,"\\tSM:", str_unique_id,"\" ", 
                                                               args.str_genome_fa, " ", args.str_sample_file_left_fq, " ", args.str_sample_file_right_fq, " > ", str_sam ] ),
                                            lstr_cur_dependencies = [ args.str_genome_fa, args.str_sample_file_left_fq, args.str_sample_file_right_fq ],
                                            lstr_cur_products = [ str_sam ] )
    
    # java -jar SortSam.jar I=Input.sam O=output.bam SO=coordinate
    cmd_bwa_bam = Command.Command( str_cur_command = "".join( [ "java -jar SortSam.jar SO=coordinate I=", str_sam, " O=", str_bam ] ),
                                            lstr_cur_dependencies = [ str_sam ],
                                            lstr_cur_products = [ str_bam ] )
    lcmd_dna_mapping_commands.extend( [ cmd_bwa_sam, cmd_bwa_bam ] )
                                      
    return { INDEX_CMD: lcmd_dna_mapping_commands, INDEX_FILE: str_bam, INDEX_FOLDER:args_call.str_file_base }


def func_do_recalibration_gatk( args_call, str_align_file, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies, logr_cur ):

    # Check for the known vcf file
    # If it does not exist, warn that the associated steps will not be ran.
    if args_call.str_vcf_file is None:
        logr_cur.warn( "".join( [ "\n\n\nWARNING, WARNING, WARNING, WARNING.\n", 
                      "GATK Recalibration: A vcf file with known variants was not provided for realignment and recalibration steps.\n",
                       "These steps may perform better if such a vcf file is provided.\n\n\n" ] ) )

    # Files
    str_dedupped_bam = os.path.join( str_tmp_dir, "dedupped.bam" )
    str_dedupped_bai = os.path.join( str_tmp_dir, "dedupped.bai" )
    str_intervals = os.path.join( str_tmp_dir, "forIndelRealigner.intervals" )
    str_qc_metrics = os.path.join( str_tmp_dir, "mark_duplicates_qc_metrics.txt" )
    str_realigned_bam = os.path.join( str_tmp_dir, "realigned.bam" )
    str_realigned_bai = os.path.join( str_tmp_dir, "realigned.bai" )
    str_recalibrated_alignment_file = os.path.join( str_tmp_dir, "recal_table.table" )
    str_recalibrated_bam = os.path.join( str_tmp_dir, "recalibrated.bam" )
    str_recalibrated_bai = os.path.join( str_tmp_dir, "recalibrated.bai" )
    str_recalibration_plots_pdf = os.path.join( str_tmp_dir, "recalibration.pdf" )
    str_sorted_bam = os.path.join( str_tmp_dir, "sorted.bam" )
    str_split_bam = os.path.join( str_tmp_dir, "split.bam" )
    str_split_bai = os.path.join( str_tmp_dir, "split.bai" )

    # This is the file that is returned, could be many of the files below depending on the settings
    # This is dynamically set and different parts of this pipeline segment are activated.
    str_return_bam = ""

    # Allows the known variants vcf file to be available or not.
    lstr_known_vcf = [] if args_call.str_vcf_file is None else [ "-known", args_call.str_vcf_file ]
    lstr_known_two_dash_vcf = [] if args_call.str_vcf_file is None else [ "--known", args_call.str_vcf_file ]

    lcmd_gatk_recalibration_commands = []
    
    # Create commands
    # SAM to BAM and QC
    cmd_add_or_replace_groups = Command.Command( str_cur_command = "".join( [ "java -jar AddOrReplaceReadGroups.jar I=", str_align_file,
                                                                     " O=", str_sorted_bam, " SO=coordinate RGID=id RGLB=library RGPL=",
                                                                     args_call.str_sequencing_platform, " RGPU=machine RGSM=", str_unique_id ] ),
                                            lstr_cur_dependencies = lstr_dependencies,
                                            lstr_cur_products = [ str_sorted_bam ] )
    cmd_mark_duplicates = Command.Command( str_cur_command = "".join( [ "java -jar MarkDuplicates.jar I=", str_sorted_bam, " O=", str_dedupped_bam, 
                                                                      " CREATE_INDEX=true M=", str_qc_metrics ] ),
                                            lstr_cur_dependencies = [ str_sorted_bam ],
                                            lstr_cur_products = [ str_dedupped_bam, str_dedupped_bai, str_qc_metrics ] )
    cmd_split_cigar_reads = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R", args_call.str_genome_fa,
                                                                      "-I", str_dedupped_bam, "-o", str_split_bam, "-rf ReassignOneMappingQuality",
                                                                      "-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS" ] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, str_dedupped_bam, str_dedupped_bai ],
                                            lstr_cur_products = [ str_split_bam, str_split_bai ] )
    lcmd_gatk_recalibration_commands.extend( [ cmd_add_or_replace_groups, cmd_mark_duplicates, cmd_split_cigar_reads ] )
    str_return_bam = str_split_bam

    # optional indel realignment step
    if not args_call.f_stop_optional_realignment:
        cmd_realigner_target_creator = Command.Command( str_cur_command = " ".join( [ "java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R",
                                                                          args_call.str_genome_fa, "-I", str_split_bam, "--out", str_intervals ] + lstr_known_two_dash_vcf ),
                                                lstr_cur_dependencies = [ args_call.str_genome_fa, str_split_bam, str_split_bai, args_call.str_vcf_file ],
                                                lstr_cur_products = [ str_intervals ] )
        cmd_indel_realigner = Command.Command( str_cur_command = " ".join( [ "java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R", 
                                                                           args_call.str_genome_fa, "-I", str_split_bam, "-targetIntervals",
                                                                           str_intervals, "--out", str_realigned_bam ] + lstr_known_vcf ),
                                                 lstr_cur_dependencies = [ args_call.str_genome_fa, str_split_bam, str_split_bai, str_intervals, args_call.str_vcf_file ],
                                                 lstr_cur_products = [ str_realigned_bam, str_realigned_bai ] )
        lcmd_gatk_recalibration_commands.extend( [ cmd_realigner_target_creator, cmd_indel_realigner ] )
        str_return_bam = str_realigned_bam

    # Recalibrate alignments
    if not args_call.str_vcf_file is None:
        cmd_base_recalibrator = Command.Command( str_cur_command = " ".join( [ "java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I", 
                                                                      str_split_bam if args_call.f_stop_optional_realignment else str_realigned_bam,
                                                                      "-R", args_call.str_genome_fa, "--out", str_recalibrated_alignment_file, "-knownSites", args_call.str_vcf_file ] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, args_call.str_vcf_file] +
                                                                  [ str_split_bam, str_split_bai ] if args_call.f_stop_optional_realignment else [ str_realigned_bam, str_realigned_bai ],
                                            lstr_cur_products = [ str_recalibrated_alignment_file ] )
        cmd_print_reads = Command.Command( str_cur_command = " ".join( [ "java -Xmx2g -jar GenomeAnalysisTK.jar -R", args_call.str_genome_fa,
                                                                       "-T PrintReads", "--out", str_recalibrated_bam, "-I", 
                                                                       str_split_bam if args_call.f_stop_optional_realignment else str_realigned_bam, 
                                                                       "--BQSR", str_recalibrated_alignment_file ] ),
                                             lstr_cur_dependencies = [ args_call.str_genome_fa, str_recalibrated_alignment_file ] +
                                                                  [ str_split_bam, str_split_bai ] if args_call.f_stop_optional_realignment else [ str_realigned_bam, str_realigned_bai ],
                                             lstr_cur_products = [ str_recalibrated_bam, str_recalibrated_bai ] )
        lcmd_gatk_recalibration_commands.extend( [ cmd_base_recalibrator, cmd_print_reads ] )
        str_return_bam = str_recalibrated_bam

        # Optional plotting of recalibration
        if args_call.f_optional_recalibration_plot:
            cmd_analyse_covariates = Command.Command( str_cur_command = " ".join( [ "java -Xmx4g -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R",
                                                                        args_call.str_genome_fa, "-BQSR", str_recalibrated_alignment_file,
                                                                        "-plots", str_recalibration_plots_pdf ] ),
                                              lstr_cur_dependencies = [ args_call.str_genome_fa, str_recalibrated_alignment_file ],
                                              lstr_cur_products = [ str_recalibration_plots_pdf ] )
            lcmd_gatk_recalibration_commands.append( cmd_analyse_covariates )

    return { INDEX_CMD : lcmd_gatk_recalibration_commands, INDEX_FILE : str_return_bam, INDEX_FOLDER : str_tmp_dir }


def func_do_rnaseq_caller_gatk( args_call, str_input_bam, str_unique_id, str_project_dir, str_tmp_dir ):
    
    # Files 
    str_filtered_variants_file = os.path.join( str_project_dir, "_".join( [ str_unique_id, "filtered_variants.vcf" ] ) )
    str_variants_file = os.path.join( str_tmp_dir, "variants.vcf" )
    str_input_bai = ".".join( [ os.path.splitext( str_input_bam )[ 0 ], "bai" ] )
    
    # Variant calling
    cmd_haplotype_caller = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R", args_call.str_genome_fa,
                                                           "-I", str_input_bam, "-recoverDanglingHeads -dontUseSoftClippedBases",
                                                           "-stand_call_conf 20.0 -stand_emit_conf 20.0 --out", str_variants_file] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, str_input_bam, str_input_bai ],
                                            lstr_cur_products = [ str_variants_file ] ).func_set_dependency_clean_level( [ str_input_bam, str_input_bai ], Command.CLEAN_NEVER )
    cmd_variant_filteration = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T VariantFiltration -R", 
                                                                     args_call.str_genome_fa, "-V", str_variants_file, "-window 35",
                                                                     "-cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD",
                                                                     "-filter \"QD < 2.0\" --out", str_filtered_variants_file ] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, str_variants_file ],
                                            lstr_cur_products = [ str_filtered_variants_file ] ).func_set_dependency_clean_level( [ str_variants_file ], Command.CLEAN_NEVER )

    return [ cmd_haplotype_caller, cmd_variant_filteration ]


def func_do_variant_calling_gatk( args_call, str_align_file, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies, logr_cur ):
    """
    Creates the commands for the GATK variant calling pipeline.
    """

    # Commands which will be returned
    lcmd_gatk_variants_commands = []

    # Perform recalibration
    dict_recalibration = func_do_recalibration_gatk( args_call, str_align_file, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies, logr_cur )
    lcmd_gatk_variants_commands.extend( dict_recalibration[ INDEX_CMD ] )
    
    # Do calling
    lcmd_rnaseq_calling = func_do_rnaseq_caller_gatk( args_call, dict_recalibration[ INDEX_FILE ], str_unique_id, str_project_dir, str_tmp_dir )    
    lcmd_gatk_variants_commands.extend( lcmd_rnaseq_calling )

    return lcmd_gatk_variants_commands


def func_call_dnaseq_like_rnaseq( args_call, str_align_file, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies, logr_cur ):
    """
    Manages the calls for calling mutations in DNA-seq data in a similar way to the RNA-seq calls.
    Development use only, not intended to be used with studies.
    Used to validate technical properties of the RNA-seq pipeline.
    
    args_call : Arguments
              : Arguments used to run pipeline
              
    Return : List of commands
           : List of commands to run for BWA alignment
    """
    
    str_dedup_bam = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "sorted.dedupped.bam" ] ) )
    str_dedup_metrics = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "sorted.dedup.metrics" ] ) )
    str_filtered_variants_file = os.path.join( str_project_dir, ".".join( [ str_unique_id, "filtered.variants.vcf" ] ) )
    str_intervals = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "religner.intervals" ] ) )
    str_raw_vcf = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "variants.vcf" ] ) )
    str_realigned_bam = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "sorted.dedup.groups.realigned.bam" ] ) )
    str_realigned_bai = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "sorted.dedup.groups.realigned.bai" ] ) )
    str_recal_plot = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "recal.pdf" ] ) )
    str_recal_snp_bam = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "recal_snp.bam" ] ) )
    str_recal_snp_bai = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "recal_snp.bai" ] ) )
    str_recal_table = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "recal.table" ] ) )
    str_recal_table_2 = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "recal_2.table" ] ) )
    str_replace_bam = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "sorted.dedup.groups.bam" ] ) )
    str_replace_bai = os.path.join( str_tmp_dir, ".".join( [ str_unique_id, "sorted.dedup.groups.bai" ] ) )


    # DNA-seq best practices
    # java -jar MarkDuplicates.jar I=input.sam O=output.bam
    cmd_dedup = Command.Command( str_cur_command = "".join( [ "java -jar MarkDuplicates.jar I=", str_align_file,
                                                             " M=", str_dedup_metrics, " O=", str_dedup_bam ] ),
                                            lstr_cur_dependencies = [ str_align_file ],
                                            lstr_cur_products = [ str_dedup_metrics, str_dedup_bam ] )
    
    # java -jar AddOrReplaceReadGroups.jar I=input.bam O=output.bam RGID=x RGLB=x RGPL=x RGPU=x RGSM=x RGCN=x RGDT=x
    cmd_replace = Command.Command( str_cur_command = "".join( [ "java -jar AddOrReplaceReadGroups.jar I=", str_dedup_bam, " O=", str_replace_bam,
                                                              " RGCN=RGCN RGID=", str_unique_id, " RGLB=library RGDT=", datetime.date.today().isoformat()," RGPL=",
                                                                args_call.str_sequencing_platform, " RGPU=machine RGSM=", str_unique_id, " CREATE_INDEX=TRUE" ] ),
                                            lstr_cur_dependencies = [ str_dedup_bam ],
                                            lstr_cur_products = [ str_replace_bam, str_replace_bai ] )
    ## Indel Realignment
    # java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R human.fasta -I original.bam -known indels.vcf -o religner.intervals
    cmd_create_target = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R",
                                                                args_call.str_genome_fa, "-I", str_replace_bam, "--out", str_intervals,
                                                                "-known", args_call.str_vcf_file ] ),
                                lstr_cur_dependencies = [ str_replace_bam, args_call.str_vcf_file ],
                                lstr_cur_products = [ str_intervals ] )
    
    # java -jar GenomeAnalysisTK.jar -T IndelRealigner -R human.fasta -I orginal.bam -known indels.vcf -targetIntervals realigner.intervals -o realigned.bam
    cmd_realign = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T IndelRealigner -R", 
                                                args_call.str_genome_fa, "-I", str_replace_bam, "-targetIntervals",
                                                str_intervals, "--out", str_realigned_bam, "-known", args_call.str_vcf_file ] ),
                                lstr_cur_dependencies = [ args_call.str_genome_fa, str_replace_bam, str_intervals ],
                                lstr_cur_products = [ str_realigned_bam, str_realigned_bai ] )
    
    ## Base Recalibration
    # java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R human.fasta -I realigned.bam -knownSites x.vcf -o recal.table
    cmd_recalibrate = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R", 
                                                args_call.str_genome_fa, "-I", str_realigned_bam, "--out", str_recal_table,
                                                "-knownSites", args_call.str_vcf_file ] ),
                                lstr_cur_dependencies = [ args_call.str_genome_fa, args_call.str_vcf_file, str_realigned_bam ],
                                lstr_cur_products = [ str_recal_table ] )
    
    # java -jar GenomeAnalysisTK.jar -T PrintReads -R human.fasta -I realigned.bam -BQSR recal.table -o recal.bam
    cmd_print = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T PrintReads -R", 
                                                args_call.str_genome_fa, "-I", str_realigned_bam, "--out", str_recal_snp_bam,
                                                "-BQSR", str_recal_table ] ),
                                lstr_cur_dependencies = [ args_call.str_genome_fa, str_realigned_bam, str_realigned_bai, str_recal_table ],
                                lstr_cur_products = [ str_recal_snp_bam, str_recal_snp_bai ] )

    ### Make plots
    # java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R human.fasta -I realigned.bam -knownSites x.vcf -BQSR recal.table -o after_recal.table
    cmd_recalibrate_2 = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R", 
                                        args_call.str_genome_fa, "-I", str_realigned_bam, "--out", str_recal_table_2,
                                        "-knownSites", args_call.str_vcf_file, "-BQSR", str_recal_table ] ),
                                lstr_cur_dependencies = [ args_call.str_genome_fa, str_realigned_bam, str_realigned_bai, 
                                                         args_call.str_vcf_file, str_recal_table ],
                                lstr_cur_products = [ str_recal_table_2 ] )

    # Commands so far
    ls_cmds = [ cmd_dedup, cmd_replace, cmd_create_target, cmd_realign, cmd_recalibrate, cmd_print, cmd_recalibrate_2 ] 

    # Optional plotting of recalibration
    if args_call.f_optional_recalibration_plot:
        # java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R human.fasta -before recal.table -after after_recal.tale -plots recal_plots.pdf
        cmd_covariates = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R", 
                                        args_call.str_genome_fa, "-before", str_recal_table, "-after", str_recal_table_2,
                                        "-plots", str_recal_plot ] ),
                                lstr_cur_dependencies = [ args_call.str_genome_fa, str_recal_table, str_recal_table_2 ],
                                lstr_cur_products = [ str_recal_plot ] )
	ls_cmds.append( cmd_covariates )
    
    # Call mutations - Single File, variant only calling in DNA-seq.
    # https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
    # java -jar GenomeAnalysisTk.jar -T HaplotypeCaller -R reference/file.fasta -I recal.bam -stand_call_conf 30 -stand_emit_conf -o output.vcf
    cmd_haplotype_caller = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R", args_call.str_genome_fa,
                                                           "-I", str_recal_snp_bam, "-stand_call_conf 30.0 -stand_emit_conf 10.0 -o", str_raw_vcf] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, str_recal_snp_bam, str_recal_snp_bai ],
                                            lstr_cur_products = [ str_raw_vcf ] ).func_set_dependency_clean_level( [ str_recal_snp_bam, str_recal_snp_bai ], Command.CLEAN_NEVER )
    
    # Hard filter like the RNA-seq
    cmd_variant_filteration = Command.Command( str_cur_command = " ".join( [ "java -jar GenomeAnalysisTK.jar -T VariantFiltration -R", 
                                                                     args_call.str_genome_fa, "-V", str_raw_vcf, "-window 35",
                                                                     "-cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD",
                                                                     "-filter \"QD < 2.0\" --out", str_filtered_variants_file ] ),
                                            lstr_cur_dependencies = [ args_call.str_genome_fa, str_raw_vcf ],
                                            lstr_cur_products = [ str_filtered_variants_file ] ).func_set_dependency_clean_level( [ str_filtered_variants_file ], Command.CLEAN_NEVER  )

    ls_cmds.extend( [ cmd_haplotype_caller, cmd_variant_filteration ] )
    return ls_cmds


def func_do_variant_calling_samtools( args_call, str_align_file, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies, logr_cur ):
    """
    Creates the commands for the SamTools variant calling pipeline.
    
    * str_reads_path : File path that points to either the Sam or Bam file
    * str_reference_file :  Fasta reference file, used in indexing the reads.
    * str_main_dir : Directory to place main results
    * str_tmp_dir : Directory to place intermediary files
    """

    # Commands to run
    lcmd_samtools_variants_commands = []
    # The bam file, stored here because the path may be changed if the optional sam -> conversion is needed.
    str_bam = str_align_file
    # Sorted bam file path
    str_bam_file = os.path.split( str_bam )[1]
    str_bam_sorted = os.path.join( str_tmp_dir, ".".join( [ os.path.splitext( str_bam_file )[0],"sorted", "bam" ] ) )
    # Index of the sorted bam file
    str_bam_sorted_index = ".".join( [ str_bam_sorted, "bai" ] )
    # Binary variant calling file
    str_bam_sorted_file = os.path.split( str_bam_sorted )[1]
    str_variants_bcf = os.path.join( str_project_dir, ".".join( [ os.path.splitext( str_bam_sorted_file )[0],"bcf" ] ) )
    # Uncompressed variant calling file
    str_variants_vcf = ".".join( [ os.path.splitext( str_bam_sorted )[0],"vcf" ] )

    # Optional SAM to BAM
    if os.path.splitext( str_align_file )[1].lower() == ".sam":
        str_bam = ".".join( [ os.path.splitext( str_align_file )[0],"bam" ] )
        lcmd_samtools_variants_commands.extend( [ 
                            Command.Command( str_cur_command = " ".join( [ "samtools view -b -S -o",str_bam, str_align_file ] ),
                                            lstr_cur_dependencies = lstr_dependencies,
                                            lstr_cur_products = [ str_bam ] ) ] )

    # Either prepare bams wit GATK best practices or minimally
    if args_call.f_recalibrate_sam:
        # Create commands for recalibration
        # Update files to recalibrated files
        dict_recalibration = func_do_recalibration_gatk( args_call, str_align_file, str_unique_id, str_project_dir, str_tmp_dir, lstr_dependencies )
        lcmd_samtools_variants_commands.extend( dict_recalibration[ INDEX_CMD ] )
        str_bam_sorted = dict_recalibration[ INDEX_FILE ]
        str_bam_sorted_index = ".".join( [ os.path.splitext( str_bam_sorted )[ 0 ], "bai" ] )
    else:
        lcmd_samtools_variants_commands.extend( [ 
                            Command.Command( str_cur_command = " ".join( [ "samtools sort", str_bam, os.path.splitext( str_bam_sorted )[0] ] ),
                                            lstr_cur_dependencies = [ str_bam ],
                                            lstr_cur_products = [ str_bam_sorted ] ),
                            Command.Command( str_cur_command = " ".join( [ "samtools index", str_bam_sorted ] ),
                                            lstr_cur_dependencies = [ str_bam_sorted ],
                                            lstr_cur_products = [ str_bam_sorted_index ] ) ] )   
    # Identify variants
    lcmd_samtools_variants_commands.extend( [ Command.Command( str_cur_command = " ".join( [ "samtools mpileup -g -f", args_call.str_genome_fa, str_bam_sorted, ">", str_variants_bcf ] ),
                                            lstr_cur_dependencies = [ str_bam_sorted ],
                                            lstr_cur_products = [ str_variants_bcf ] ),
                            Command.Command( str_cur_command = " ".join( [ "bcftools view -vcg", str_variants_bcf, ">", str_variants_vcf ] ),
                                            lstr_cur_dependencies = [ str_variants_bcf ],
                                            lstr_cur_products = [ str_variants_vcf ] )  ] )
    return lcmd_samtools_variants_commands


def run( args_call, f_do_index = False ):
    """
    Runs the pipeline. This is placed in a function so that multiple
    scripts with different arguments requirements can run it.
    For instance, one script managing the complete pipeline requiring many more arguments
    while having a simple script running only the first indexing step to make a global index
    useable in all subsequent runs associated with the reference genome generating the index.
    
    * args_call : Arguments
                : Arguments used to run pipeline.
                
    * f_do_index : Boolean
                 : True indicates ONLY the index will be performed.
    """

    # Reset args if validating with DNA-seq data
    if args.f_validate_by_dnaseq:
        args.str_alignment_mode = STR_DNASEQ_VALIDATION
        args.str_variant_call_mode = STR_DNASEQ_VALIDATION
        args.f_recalibrate_sam = True

    # Constants
    # If not using a premade index and indexing only, do not update the name of the index dir with the sample name
    # Does not need to be unique. Otherwise update the name of the index directory so it is unique,
    # In case the pipeline is ran for multiple samples at once.
    str_sample_postfix = os.path.splitext( os.path.basename( args_call.str_genome_fa if f_do_index 
                                                             else args_call.str_sample_file_left_fq ) )[ 0 ]
    str_sample_postfix = str_sample_postfix.replace(".","_")

    STR_MISC_DIR = "_".join( [ "misc", str_sample_postfix ] )

    # Make the gtf files required for TOPHAT and GSNAP flows
    if args_call.str_alignment_mode in [ STR_ALIGN_GSNP, STR_ALIGN_TOPHAT ]:
        if not args_call.str_gtf_file_path:
            print "".join( [ "GTF file is required when using the ",args_call.str_alignment_mode,
                             " alignment method. Please provide and try again." ] )
            exit( 5 )

    # Vary alignment depending on arguments
    dict_align_funcs = { STR_ALIGN_GSNP : func_do_gsnp_alignment,
                        STR_ALIGN_STAR_LIMITED : func_do_star_alignment,
                        STR_ALIGN_STAR : func_do_star_alignment,
                        STR_ALIGN_TOPHAT : func_do_top_hat_alignment,
                        STR_DNASEQ_VALIDATION : func_do_BWA_alignment }
    
    # Vary variant calling depending on arguments
    dict_variant_calling_funcs = { STR_DNASEQ_VALIDATION : func_call_dnaseq_like_rnaseq,
                                   STR_VARIANT_GATK : func_do_variant_calling_gatk,
                                   STR_VARIANT_SAMTOOLS : func_do_variant_calling_samtools }

    # If the output directory is not given, get the file base from a sample file
    if not args_call.str_file_base:
        args_call.str_file_base = os.path.splitext( os.path.basename( args_call.str_sample_file_left_fq ) )[0]
        args_call.str_file_base = args_call.str_file_base.replace(".","_")

    # Make sure the output directory is absolute
    args_call.str_file_base = os.path.abspath( args_call.str_file_base )

    # Base outputs on the sample file unless an output directory is given
    # Directories
    str_misc_dir = os.path.join( args_call.str_file_base, STR_MISC_DIR )

    # Make pipeline object and indicate Log file
    pline_cur = Pipeline.Pipeline( str_name = "rnaseq_mutation", 
                                   str_log_to_file = args_call.str_log_file, 
                                   str_update_source_path = args_call.str_update_classpath if hasattr( args_call, "str_update_classpath" ) else None )

    # Put pipeline in test mode if needed.
    if args_call.f_Test:
        pline_cur.func_test_mode()

    # Make commands bsub if indicated
    if args_call.str_bsub_queue:
        pline_cur.func_do_bsub( str_memory = args_call.str_max_memory, str_queue = args_call.str_bsub_queue )
    
    # Check for installed software
    #if not pline_cur.func_check_installed( [ "STAR" ]  ):
    #    exit( 1 )
    
    # Check to make sure input files exist
    lstr_files_to_check = [ args_call.str_genome_fa ]
    if not f_do_index:
        lstr_files_to_check.extend( Command.Command.func_remove_temp_files( [ args_call.str_sample_file_left_fq, args_call.str_sample_file_right_fq ] ))
    if not pline_cur.func_check_files_exist( lstr_files_to_check ):
        if not args_call.f_Test:
            exit( 2 )

    # Handle indexing and alignment
    # Vary handling based on alignment type
    dict_align_info = dict_align_funcs[ args_call.str_alignment_mode ]( args_call = args_call,
                                                                        str_unique_id = str_sample_postfix, 
                                                                        pline_cur = pline_cur, f_index_only = f_do_index )
    # Commands to be ran, here the alignment commands
    lcmd_commands = dict_align_info[ INDEX_CMD ]
        
    # Just do the initial indexing. No other part
    if f_do_index:
        # Run commands lcmd_commands, str_output_dir, i_clean_level = Command.CLEAN_NEVER, str_run_name = ""
        if not pline_cur.func_run_commands( lcmd_commands = lcmd_commands, str_output_dir = args_call.str_file_base, f_clean = args_call.f_clean ):
            exit( 4 )
        exit( 0 )
        
    # Make directories
    if not pline_cur.func_mkdirs( [ str_misc_dir ] ):
        exit( 3 )

    # Add variant calling commands
    lcmd_commands.extend( dict_variant_calling_funcs[ args_call.str_variant_call_mode ]( args_call = args_call,
                                                                                        str_align_file = dict_align_info[ INDEX_FILE ],
                                                                                        str_unique_id = str_sample_postfix,
                                                                                        str_project_dir = args_call.str_file_base, 
                                                                                        str_tmp_dir = str_misc_dir, 
                                                                                        lstr_dependencies = [ dict_align_info[ INDEX_FOLDER ] ],
                                                                                        logr_cur = pline_cur.logr_logger ) )

    # Run commands
    if not pline_cur.func_run_commands( lcmd_commands = lcmd_commands, str_output_dir = args_call.str_file_base, f_clean = args_call.f_clean ):
        exit( 4 )
    exit( 0 )


if __name__ == "__main__":

    import argparse
        
    # Parse arguments
    prsr_arguments = argparse.ArgumentParser( prog = "rnaseq_mutation_pipeline.py", description = "Variant calling using RNASeq NGS sequencing", formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    prsr_arguments.add_argument( "-a", "--realign", dest = "f_stop_optional_realignment", default = False, action = "store_true", help = "Turns off optional indel realignment step." )
    prsr_arguments.add_argument( "-b", "--bsub_queue", metavar = "BSUB_Queue", dest = "str_bsub_queue", default = None, help = "If given, each command will sequentially be ran on this queue with bsub." )
    prsr_arguments.add_argument( "-c", "--clean", dest = "f_clean", default = False, action="store_true", help = "Turns on (true) or off (false) cleaning of intermediary product files." ) 
    prsr_arguments.add_argument( "-d", "--alignment_mode", metavar = "Alignment_mode", dest = "str_alignment_mode", default = STR_ALIGN_STAR, choices = LSTR_ALIGN_CHOICES, help = "Specifies the alignment and indexing algorithm to use." )
    prsr_arguments.add_argument( "-e", "--variant_call_mode", metavar = "Call_mode", dest = "str_variant_call_mode", default = STR_VARIANT_GATK, choices = LSTR_VARIANT_CALLING_CHOICES, help = "Specifies the variant calling method to use." )
    prsr_arguments.add_argument( "-f", "--reference", metavar = "Reference_genome", dest = "str_genome_fa", required = True, help = "Path to the reference genome to use in the analysis pipeline." )
    prsr_arguments.add_argument( "-g", "--log", metavar = "Optional_logging_file", dest = "str_log_file", default = None, help = "Optional log file, if not given logging will be to the standard out." )
    prsr_arguments.add_argument( "-i", "--index", metavar = "Use_premade_index", dest = "str_initial_index", default = None, help = "The initial index is made only from the reference genome and can be shared. If premade, supply a path here to the index directory so that it is not rebuilt for every alignment. Please provide the full path." )
    prsr_arguments.add_argument( "-j", "--recalibrate_sam", dest = "f_recalibrate_sam", default = True, action="store_false", help = "If used, turns off gatk recalibration of bam files before samtools variant calling." ) 
    prsr_arguments.add_argument( "-k", "--gtf", metavar = "Reference GTF", dest = "str_gtf_file_path", default = None, help = "GTF file for reference genome.")
    prsr_arguments.add_argument( "-l", "--left", metavar = "Left_sample_file", dest = "str_sample_file_left_fq", required = True, help = "Path to one of the two paired RNAseq samples ( left )" )
    prsr_arguments.add_argument( "-m", "--max_bsub_memory", metavar = "Max_BSUB_Mem", dest = "str_max_memory", default = "8", help = "The max amount of memory in GB requested when running bsub commands." )
    prsr_arguments.add_argument( "-n", "--threads", metavar = "Process_threads", dest = "i_number_threads", type = int, default = 1, help = "The number of threads to use for multi-threaded steps." )
    prsr_arguments.add_argument( "-o", "--out_dir", metavar = "Output_directory", dest = "str_file_base", default = None, help = "The output directory where results will be placed. If not given a directory will be created from sample names and placed with the samples." )
    prsr_arguments.add_argument( "-p", "--plot", dest = "f_optional_recalibration_plot", default = True, action = "store_false", help = "Turns off plotting recalibration of alignments." )
    prsr_arguments.add_argument( "-r", "--right", metavar = "Right_sample_file", dest = "str_sample_file_right_fq", required = True, help = "Path to one of the two paired RNAseq samples ( right )" )
    prsr_arguments.add_argument( "-s", "--sequencing_platform", metavar = "Sequencing Platform", dest = "str_sequencing_platform", default = "ILLUMINA", choices = LSTR_SEQ_CHOICES, help = "The sequencing platform used to generate the samples choices include " + " ".join( LSTR_SEQ_CHOICES ) + "." )
    prsr_arguments.add_argument( "-t", "--test", dest = "f_Test", default = False, action = "store_true", help = "Will check the environment and display commands line but not run.")
    prsr_arguments.add_argument( "-u", "--update_command", dest = "str_update_classpath", default = None, help = "Allows a class path to be added to the jars. eg. 'command.jar:/APPEND/THIS/PATH/To/JAR,java.jar:/Append/Path'")
    prsr_arguments.add_argument( "--validate_dnaseq", dest = "f_validate_by_dnaseq", default = False, action = "store_true", help = "Used for development only. Should not be used with biological samples.")
    prsr_arguments.add_argument( "-w", "--vcf", metavar = "Variant_calling_file_for_the_reference_genome", dest = "str_vcf_file", default = None, help = "Variant calling file for the reference genome.")
    prsr_arguments.add_argument( "-y", "--star_memory", metavar = "Star_memory", dest = "str_star_memory_limit", default = None, help = "Memory limit for star index. This should be used to increase memory if needed. Reducing memory consumption should be performed with the STAR Limited mod." )
    args = prsr_arguments.parse_args()
    
    run( args )
