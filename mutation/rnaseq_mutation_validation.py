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
STR_RNA_BAM = "RNA_BAM"
STR_RNA_SAMPLE_DEPTH = "RNA_DEPTH"
STR_RNA_VCF_SNP = "RNA_VCF_SNP"
STR_DNA_BAM = "DNA_BAM"
STR_DNA_SAMPLE_DEPTH = "DNA_DEPTH"
STR_DNA_VCF = "DNA_VCF"
STR_DNA_VCF_SNP = "DNA_VCF_SNP"
STR_MAF_SAMPLE_NAME = "MAF_SAMPLE"

# Paired sample file indices
I_RNA_VCF = 4
I_RNA_BAM = 3
I_DNA_VCF = 2
I_DNA_BAM = 1
I_MAF_SAMPLE = 0

class RNASEQ_mutation_validation( ParentScript.ParentScript ):

    def __init__( self ):

        # File structure
        self.str_depth_dir = "depth"
        self.str_figure_dir = "figure"
        self.str_igv = "igv"
        self.str_log_dir = "log"
        self.str_src_dir = "src"
        self.str_tab_dir = "tab"
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

        dict_return = {}
        f_skip_header = True
        with open( str_paired_file, "r" ) as hndl_paired_file:
          for lstr_line in csv.reader( hndl_paired_file, delimiter = "\t" ):
            if f_skip_header:
              f_skip_header = False
              continue
            dict_return[ lstr_line[ I_RNA_VCF ]] = { STR_DNA_BAM:lstr_line[ I_DNA_BAM ],
                                                   STR_RNA_BAM:lstr_line[ I_RNA_BAM ],
                                                   STR_DNA_VCF:lstr_line[ I_DNA_VCF ],
                                                   STR_MAF_SAMPLE_NAME:lstr_line[ I_MAF_SAMPLE ] }
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

        # File that pairs RNA and DNA seq files
        # dest = str_sample_pairing_file
        arg_raw.add_argument( "-s", "--sample_file", dest = "str_sample_pairing_file", required = True, help = "A file that contains pairing of RNA and DNA Bam and Vcf files." )

        # MAF file
        # dest = str_maf_file
        arg_raw.add_argument( "-m", "--maf", dest = "str_maf_file", default=None, help = "The maf file used to compare DNA dn RNA vcfs." )

        # Reduced dbSNP vcf file
        # dest = str_reduced_dbsnp_vcf
        arg_raw.add_argument( "-d" , "--ref_vcf", dest = "str_reduced_dbsnp_vcf", required = True, help = "A reference vcf file." )

        # Expression matrix
        # dest = str_expression_matrix
        arg_raw.add_argument( "-x", "--expression_matrix", dest = "str_expression_matrix", required = True, help = "Expression matrix from the RNA samples." )

        # GTF file for the current reference annotation
        # dest = str_annotation_gtf
        arg_raw.add_argument( "-g", "--gtf", dest = "str_annotation_gtf", required = True, help = "GTF file matching the referenc RNA samples." )

        # Comma delimited string of key mutations to look at occurence in samples
        # dest = args_parsed.str_key_mutations
        arg_raw.add_argument( "--key", dest = "str_key_mutations", required = True, help = "Comma delimited string of key mutations (gene names) to look at occurence in samples." )

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """

        * return : List of commands to be ran
                 : list of command objects
        """
        # Check to make sure directories exit
        STR_DEPTH_DIR = os.path.join( args_parsed.str_file_base, self.str_depth_dir )
        STR_FIGURE_DIR = os.path.join( args_parsed.str_file_base, self.str_figure_dir )
        STR_LOG_DIR = os.path.join( args_parsed.str_file_base, self.str_log_dir )
        STR_TAB_DIR = os.path.join( args_parsed.str_file_base, self.str_tab_dir )
        STR_DNA_FIGURE = os.path.join( args_parsed.str_file_base, self.str_figure_dir, self.str_figure_dna )
        STR_RNA_FIGURE = os.path.join( args_parsed.str_file_base, self.str_figure_dir, self.str_figure_rna )
        STR_MAF_FIGURE = os.path.join( args_parsed.str_file_base, self.str_figure_dir, self.str_figure_maf )
        STR_FILTERED_VCFS = os.path.join( args_parsed.str_file_base, self.str_filtered_vcf )
        # IGV files
        str_IGV_dir = os.path.join( STR_FIGURE_DIR, self.str_igv )
        STR_IGV_MAF_DNA = os.path.join( str_IGV_dir, "MAF_DNA" )
        STR_IGV_DNA_NOT_RNA = os.path.join( str_IGV_dir, "DNA_NOT_RNA" )
        STR_IGV_MAF_NOT_DNA = os.path.join( str_IGV_dir, "MAF_NOT_DNA" )
        STR_IGV_MAF_RNA = os.path.join( str_IGV_dir, "MAF_RNA" )
        STR_IGV_MAF_NOT_RNA = os.path.join( str_IGV_dir, "MAF_NOT_RNA" )
        STR_IGV_DNA_RNA_NOT_MAF = os.path.join( str_IGV_dir, "DNA_RNA_NOT_MAF" )
        STR_IGV_MAF_RNA_NOT_DNA = os.path.join( str_IGV_dir, "MAF_RNA_NOT_DNA" )
        STR_IGV_MAF_DNA_RNA = os.path.join( str_IGV_dir, "MAF_DNA_RNA" )
        STR_IGV_ERROR = os.path.join( str_IGV_dir, "ERROR" )
        for str_dir_path in [ STR_FIGURE_DIR, STR_LOG_DIR, STR_TAB_DIR, STR_DNA_FIGURE, STR_RNA_FIGURE, STR_MAF_FIGURE, STR_DEPTH_DIR,
                              STR_FILTERED_VCFS, str_IGV_dir, STR_IGV_MAF_DNA, STR_IGV_MAF_NOT_DNA, STR_IGV_MAF_RNA, STR_IGV_MAF_NOT_RNA,
                              STR_IGV_DNA_RNA_NOT_MAF, STR_IGV_MAF_RNA_NOT_DNA, STR_IGV_MAF_DNA_RNA, STR_IGV_ERROR, STR_IGV_DNA_NOT_RNA ]:
          if not os.path.exists( str_dir_path ):
            os.mkdir( str_dir_path )

        # Get the sample pairing file for the study
        dict_sample_study = self.parse_paired_sample_file( str_paired_file = args_parsed.str_sample_pairing_file )
        lcmd_commands = []
        lstr_maf_dna_tab = []
        lstr_maf_rna_tab = []
        lstr_dna_rna_tab = []
        lstr_rna_vcfs_snps = []
        lstr_dna_vcfs_snps = []

        # Distance metric functions
        str_jaccard_distance = os.path.join( self.str_src_dir, "jaccard_distance.R" )
        str_jaccard_distance_restricted = os.path.join( self.str_src_dir, "jaccard_distance_restricted.R" )

        # SOFT Filter the DBSNP vcf file to SNP
        str_filtered_dbsnp_vcf = os.path.join( STR_FILTERED_VCFS, os.path.splitext( os.path.basename( args_parsed.str_reduced_dbsnp_vcf ) )[ 0 ] + "_snp.vcf" )
        str_filtered_dbsnp_vcf_command = os.path.join( self.str_src_dir, "reduce_vcf_to_snps.py --reference " ) + args_parsed.str_reduced_dbsnp_vcf + " " + str_filtered_dbsnp_vcf
        lcmd_commands.append( Command.Command( str_cur_command = str_filtered_dbsnp_vcf_command,
                                               lstr_cur_dependencies = [ args_parsed.str_reduced_dbsnp_vcf ],
                                               lstr_cur_products = [ str_filtered_dbsnp_vcf ] ) )

        # Create secondary files associate with each sample (vcf) including depth files and tab files.
        for str_sample_vcf in dict_sample_study:

            dict_sample_pairing = dict_sample_study[ str_sample_vcf ]
            str_sample_key = os.path.splitext( os.path.basename( str_sample_vcf ) )[ 0 ]

            # Create depth file RNA
            str_RNA_depth_file = os.path.splitext( os.path.basename( dict_sample_pairing[ STR_RNA_BAM ] ) )[ 0 ] + "_bam.depth"
            str_RNA_depth_file = os.path.join( STR_DEPTH_DIR, str_RNA_depth_file )
            lcmd_commands.append( Command.Command( str_cur_command = "samtools depth " + dict_sample_pairing[ STR_RNA_BAM ]+ " > " + str_RNA_depth_file,
                                                   lstr_cur_dependencies = [ dict_sample_pairing[ STR_RNA_BAM ] ],
                                                   lstr_cur_products = [ str_RNA_depth_file ] ) )

            # Compress depth file RNA
            str_RNA_depth_compressed_file = os.path.splitext( str_RNA_depth_file )[ 0 ] + ".depth.gz"
            lcmd_commands.append( Command.Command( str_cur_command = "gzip " + str_RNA_depth_file,
                                                   lstr_cur_dependencies = [ str_RNA_depth_file ],
                                                   lstr_cur_products = [ str_RNA_depth_compressed_file ] ) )

            # Create depth file DNA
            str_DNA_depth_file = os.path.splitext( os.path.basename( dict_sample_pairing[ STR_DNA_BAM ] ) )[ 0 ] + "_bam.depth"
            str_DNA_depth_file = os.path.join( STR_DEPTH_DIR, str_DNA_depth_file )
            lcmd_commands.append( Command.Command( str_cur_command = "samtools depth " + dict_sample_pairing[ STR_DNA_BAM ]+ " > " + str_DNA_depth_file,
                                                   lstr_cur_dependencies = [ dict_sample_pairing[ STR_DNA_BAM ] ],
                                                   lstr_cur_products = [ str_DNA_depth_file ] ) )

            # Compress depth file DNA
            str_DNA_depth_compressed_file = os.path.splitext( str_DNA_depth_file )[ 0 ] + ".depth.gz"
            lcmd_commands.append( Command.Command( str_cur_command = "gzip " + str_DNA_depth_file,
                                                   lstr_cur_dependencies = [ str_DNA_depth_file ],
                                                   lstr_cur_products = [ str_DNA_depth_compressed_file ] ) )

            # Filter RNA and DNA VCF files to SNPS
            # SOFT Filter DNA VCF down to snps
            str_filtered_dna_vcf = os.path.join( STR_FILTERED_VCFS, os.path.basename( os.path.splitext( dict_sample_pairing[ STR_DNA_VCF ] )[ 0 ] + "_snp.vcf" ) )
            str_filtered_dna_vcf_command = os.path.join( self.str_src_dir, "reduce_vcf_to_snps.py " ) + dict_sample_pairing[ STR_DNA_VCF ] + " " + str_filtered_dna_vcf
            lcmd_commands.append( Command.Command( str_cur_command = str_filtered_dna_vcf_command,
                                                   lstr_cur_dependencies = [ dict_sample_pairing[ STR_DNA_VCF ] ],
                                                   lstr_cur_products = [ str_filtered_dna_vcf ] ) )

            # Store filtered DNA vcf
            dict_sample_pairing[ STR_DNA_VCF_SNP ] = str_filtered_dna_vcf

            # SOFT Filter RNA VCF down to snps
            str_filtered_rna_vcf = os.path.join( STR_FILTERED_VCFS, os.path.basename( os.path.splitext( str_sample_vcf )[ 0 ] + "_snp.vcf" ) )
            str_filtered_rna_vcf_command = os.path.join( self.str_src_dir, "reduce_vcf_to_snps.py " ) + str_sample_vcf + " " + str_filtered_rna_vcf
            lcmd_commands.append( Command.Command( str_cur_command = str_filtered_rna_vcf_command,
                                                   lstr_cur_dependencies = [ str_sample_vcf ],
                                                   lstr_cur_products = [ str_filtered_rna_vcf ] ) )

            # Store filtered RNA vcf
            dict_sample_pairing[ STR_RNA_VCF_SNP ] = str_filtered_rna_vcf

            # Create tab files
            # Tab files for vcf files that include if common or somatic
            str_file_maf_dna_tab = os.path.join( STR_TAB_DIR, str_sample_key + "_maf_dna.tab" )
            str_log_maf_dna = os.path.join( STR_LOG_DIR, str_sample_key + "_maf_dna.log" )
            if not args_parsed.str_maf_file is None:
                str_cmd_tab_MAF_DNA = os.path.join( self.str_src_dir, "vcfs_to_snp_calls_tab.py" ) + " --maf_reference " + args_parsed.str_maf_file + " --tumor " + dict_sample_pairing[ STR_MAF_SAMPLE_NAME ] + " --count_reference " + str_DNA_depth_compressed_file + " --vcf " + dict_sample_pairing[ STR_DNA_VCF_SNP ] + " --count " + str_DNA_depth_compressed_file + " " + str_file_maf_dna_tab + " > " + str_log_maf_dna
                lcmd_commands.append( Command.Command( str_cur_command = str_cmd_tab_MAF_DNA,
                                                   lstr_cur_dependencies = [ args_parsed.str_maf_file, str_DNA_depth_compressed_file, dict_sample_pairing[ STR_DNA_VCF_SNP ] ],
                                                   lstr_cur_products = [ str_file_maf_dna_tab, str_log_maf_dna ] ) )

                str_file_maf_rna_tab = os.path.join( STR_TAB_DIR, str_sample_key + "_maf_rna.tab" )
                str_log_maf_rna = os.path.join( STR_LOG_DIR, str_sample_key + "_maf_rna.log" )
                str_cmd_tab_MAF_RNA = os.path.join( self.str_src_dir, "vcfs_to_snp_calls_tab.py" ) + " --maf_reference " + args_parsed.str_maf_file + " --tumor " + dict_sample_pairing[ STR_MAF_SAMPLE_NAME ] + " --count_reference " + str_DNA_depth_compressed_file + " --vcf " + dict_sample_pairing[ STR_RNA_VCF_SNP ] + " --count " + str_RNA_depth_compressed_file + " " + str_file_maf_rna_tab + " > " + str_log_maf_rna
                lcmd_commands.append( Command.Command( str_cur_command = str_cmd_tab_MAF_RNA,
                                                   lstr_cur_dependencies = [ args_parsed.str_maf_file, str_DNA_depth_compressed_file, str_RNA_depth_compressed_file, dict_sample_pairing[ STR_RNA_VCF_SNP ] ],
                                                   lstr_cur_products = [ str_file_maf_rna_tab, str_log_maf_rna ] ) )

                # Store tab file for later
                lstr_maf_dna_tab.append( str_file_maf_dna_tab )
                lstr_maf_rna_tab.append( str_file_maf_rna_tab )

            str_file_dna_rna_tab = os.path.join( STR_TAB_DIR, str_sample_key + "_dna_rna.tab" )
            str_log_dna_rna = os.path.join( STR_LOG_DIR, str_sample_key + "_dna_rna.log" )
            str_cmd_tab_DNA_RNA = os.path.join( self.str_src_dir, "vcfs_to_snp_calls_tab.py" ) + " --vcf_reference " + dict_sample_pairing[ STR_DNA_VCF_SNP ] + " --count_reference " + str_DNA_depth_compressed_file + " --vcf " + dict_sample_pairing[ STR_RNA_VCF_SNP ] + " --count " + str_RNA_depth_compressed_file + " " + str_file_dna_rna_tab + " > " + str_log_dna_rna
            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_tab_DNA_RNA,
                                                   lstr_cur_dependencies = [ dict_sample_pairing[ STR_DNA_VCF_SNP ], str_DNA_depth_compressed_file, str_RNA_depth_compressed_file, dict_sample_pairing[ STR_RNA_VCF_SNP ] ],
                                                   lstr_cur_products = [ str_file_dna_rna_tab, str_log_dna_rna ] ) )

            # Store tab file for later
            lstr_dna_rna_tab.append( str_file_dna_rna_tab )

            # Store vcf files
            lstr_rna_vcfs_snps.append( dict_sample_pairing[ STR_RNA_VCF_SNP ] )
            lstr_dna_vcfs_snps.append( dict_sample_pairing[ STR_DNA_VCF_SNP ] )

            # Plot the rate of mutation over a window.
            # Mutation rates are often local.
            # Would be interesting to filter on # of mutations in a gene given the mutation rates locally.

            # Make IGV figures for interesting features
#            str_maf_dna_out = os.path.join( STR_IGV_MAF_DNA, "igv_maf_dna.txt" )
#            str_dna_not_rna_out = os.path.join( STR_IGV_DNA_NOT_RNA, "igv_dna_not_rna.txt" )
#            str_maf_not_dna_out = os.path.join( STR_IGV_MAF_NOT_DNA, "igv_maf_not_dna.txt" )
#            str_maf_rna_out = os.path.join( STR_IGV_MAF_RNA, "igv_maf_rna.txt" )
#            str_maf_not_rna_out = os.path.join( STR_IGV_MAF_NOT_RNA, "igv_maf_not_rna.txt" )
#            str_dna_rna_not_maf_out = os.path.join( STR_IGV_DNA_RNA_NOT_MAF, "igv_dna_rna_not_maf.txt" )
#            str_maf_rna_not_dna_out = os.path.join( STR_IGV_MAF_RNA_NOT_DNA, "igv_maf_rna_not_dna.txt" )
#            str_maf_dna_rna_out = os.path.join( STR_IGV_MAF_DNA_RNA, "igv_maf_dna_rna.txt" )
#            str_error_depth_out = os.path.join( STR_IGV_ERROR, "igv_depth_error.txt" )
#            str_error_out = os.path.join( STR_IGV_ERROR, "igv_error.txt" )
#            str_cmd_IGV = " ".join( [ os.path.join( self.str_src_dir, "tabs_comparisons_to_igv_batch.py" ),
#                                      "--maf_dna " + str_file_maf_dna_tab,
#                                      "--maf_rna " + str_file_maf_rna_tab,
#                                      "--dna_rna " + str_file_dna_rna_tab,
#                                      "--maf_bam " + dict_sample_pairing[ STR_DNA_BAM ],
#                                      "--dna_bam " + dict_sample_pairing[ STR_DNA_BAM ],
#                                      "--rna_bam " + dict_sample_pairing[ STR_RNA_BAM ],
#                                      "--maf_dna_out " + str_maf_dna_out,
#                                      "--dna_not_rna_out " + str_dna_not_rna_out,
#                                      "--maf_not_dna_out " + str_maf_not_dna_out,
#                                      "--maf_rna_out " + str_maf_rna_out,
#                                      "--maf_not_rna_out " + str_maf_not_rna_out,
#                                      "--dna_rna_not_maf_out " + str_dna_rna_not_maf_out,
#                                      "--maf_rna_not_dna_out " + str_maf_rna_not_dna_out,
#                                      "--maf_dna_rna_out " + str_maf_dna_rna_out,
#                                      "--error_depth_out " + str_error_depth_out,
#                                      "--error_out " + str_error_out ] )
#            lcmd_commands.append( Command.Command( str_cur_command = str_cmd_IGV,
#                                                   lstr_cur_dependencies = [ str_file_maf_dna_tab, str_file_maf_rna_tab, str_file_dna_rna_tab,
#                                                   dict_sample_pairing[ STR_DNA_BAM ], dict_sample_pairing[ STR_DNA_BAM ], dict_sample_pairing[ STR_RNA_BAM ] ],
#                                                   lstr_cur_products = [ str_maf_dna_out, str_maf_not_dna_out, str_maf_rna_out,
#                                                                         str_maf_not_rna_out, str_dna_rna_not_maf_out, str_maf_rna_not_dna_out,
#                                                                         str_maf_dna_rna_out, str_dna_not_rna_out ] ) )

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
