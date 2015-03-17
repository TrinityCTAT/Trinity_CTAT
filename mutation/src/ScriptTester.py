
__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2015"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@broadinstitute.org"
__status__ = "Development"

import Commandline
import os
import ParentPipelineTester
import unittest

class ScriptTester( ParentPipelineTester.ParentPipelineTester ):
    """
    Testing for scripts, starting at command line.
    """
    
    str_script_dir = "src"
    str_test_data_dir = "test_data"
    str_test_data_dir_working = os.path.join( str_test_data_dir, "active_testing_script_tester" )

    def not_test_reduce_vcf_to_snp_for_small_file_no_filter( self ):
        """
        Test reducing the vcf file to snps for a file that is small (less than 100) and should not be filtered.
        """

        # Create test environment
        str_filtered_vcf_script = os.path.join( self.str_script_dir, "reduce_vcf_to_snps.py" )
        str_filtered_vcf_test_file = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_small_file_no_filter.vcf" )
        str_filtered_vcf = os.path.join( self.str_test_data_dir_working, "test_reduce_vcf_to_snp_for_small_file_no_filter_RESULT.vcf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )

        # Call Example script
        str_command = " ".join( [ str_filtered_vcf_script, "--reference", str_filtered_vcf_test_file, str_filtered_vcf ] )
        Commandline.Commandline().func_CMD( str_command )
        
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_filtered_vcf, str_filtered_vcf_test_file )

        # Destroy environment
        self.func_remove_files( [ str_filtered_vcf ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )

        # Evaluate
        self.func_test_true( f_success )

    def not_test_reduce_vcf_to_snp_for_large_file_no_filter( self ):
        """
        Test reducing the vcf file to snps for a file that is large (3003) and should not be filtered.
        This is making sure the buffering in the script is working.
        """

        # Create test environment
        str_filtered_vcf_script = os.path.join( self.str_script_dir, "reduce_vcf_to_snps.py" )
        str_filtered_vcf_test_file = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_large_file_no_filter.vcf" )
        str_filtered_vcf = os.path.join( self.str_test_data_dir_working, "test_reduce_vcf_to_snp_for_large_file_no_filter_RESULT.vcf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )

        # Call Example script
        str_command = " ".join( [ str_filtered_vcf_script, "--reference", str_filtered_vcf_test_file, str_filtered_vcf ] )
        Commandline.Commandline().func_CMD( str_command )
        
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_filtered_vcf, str_filtered_vcf_test_file )

        # Destroy environment
        self.func_remove_files( [ str_filtered_vcf ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )

        # Evaluate
        self.func_test_true( f_success )

    def not_test_reduce_vcf_to_snp_for_small_file_filter_reference( self ):
        """
        Test reducing the vcf file to snps for a file that is small (less than 100) and should be filtered.
        When filtering as a reference, it will not have PASS info so this is ignored.
        """

        # Create test environment
        str_filtered_vcf_script = os.path.join( self.str_script_dir, "reduce_vcf_to_snps.py" )
        str_filtered_vcf_test_file = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_small_file_filter.vcf" )
        str_filtered_vcf_answer = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_small_file_filter_reference_ANSWER.vcf" )
        str_filtered_vcf_result = os.path.join( self.str_test_data_dir_working, "test_reduce_vcf_to_snp_for_small_file_filter_reference_RESULT.vcf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )

        # Call Example script
        str_command = " ".join( [ str_filtered_vcf_script, "--reference", str_filtered_vcf_test_file, str_filtered_vcf_result ] )
        Commandline.Commandline().func_CMD( str_command )

        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_filtered_vcf_answer, str_filtered_vcf_result )

        # Destroy environment
        self.func_remove_files( [ str_filtered_vcf_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )

        # Evaluate
        self.func_test_true( f_success )

    def not_test_reduce_vcf_to_snp_for_small_file_filter( self ):
        """
        Test reducing the vcf file to snps for a file that is small (less than 100) and should be filtered.
        """

        # Create test environment
        str_filtered_vcf_script = os.path.join( self.str_script_dir, "reduce_vcf_to_snps.py" )
        str_filtered_vcf_test_file = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_small_file_filter.vcf" )
        str_filtered_vcf_answer = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_small_file_filter_ANSWER.vcf" )
        str_filtered_vcf_result = os.path.join( self.str_test_data_dir_working, "test_reduce_vcf_to_snp_for_small_file_filter_RESULT.vcf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )

        # Call Example script
        str_command = " ".join( [ str_filtered_vcf_script, str_filtered_vcf_test_file, str_filtered_vcf_result ] )
        Commandline.Commandline().func_CMD( str_command )
        
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_filtered_vcf_answer, str_filtered_vcf_result )

        # Destroy environment
        self.func_remove_files( [ str_filtered_vcf_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )

        # Evaluate
        self.func_test_true( f_success )

    def not_test_reduce_vcf_to_snp_for_large_file_filter( self ):
        """
        Test reducing the vcf file to snps for a file that is large (less than 5000) and should be filtered.
        """

        # Create test environment
        str_filtered_vcf_script = os.path.join( self.str_script_dir, "reduce_vcf_to_snps.py" )
        str_filtered_vcf_test_file = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_large_file_filter.vcf" )
        str_filtered_vcf_answer = os.path.join( self.str_test_data_dir, "test_reduce_vcf_to_snp_for_large_file_filter_ANSWER.vcf" )
        str_filtered_vcf_result = os.path.join( self.str_test_data_dir_working, "test_reduce_vcf_to_snp_for_large_file_filter_RESULT.vcf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )

        # Call Example script
        str_command = " ".join( [ str_filtered_vcf_script, str_filtered_vcf_test_file, str_filtered_vcf_result ] )
        Commandline.Commandline().func_CMD( str_command )
        
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_filtered_vcf_answer, str_filtered_vcf_result )

        # Destroy environment
        self.func_remove_files( [ str_filtered_vcf_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )

        # Evaluate
        self.func_test_true( f_success )

    def not_test_tabs_to_percent_mutations_for_one_tab_file( self ):
        """
        Test tabs_to_percent_mutations.py for reading in one tab file with mutations.
        """
        # Create test environment
        str_count_mutations_script = os.path.join( self.str_script_dir, "tabs_to_percent_mutations.py" )
        str_count_mutations_gtf_file = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations.gtf" )
        str_count_mutations_key_genes = "Gene1,Gene3,Gene5"
        str_count_mutations_input_tab_file = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_1.tab" )
        str_count_mutations_output_file_answer = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_for_one_tab_file_ANSWER.txt" )
        str_count_mutations_output_file_result = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_one_tab_file_RESULT.txt" )
        str_count_mutations_output_file_pdf = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_one_tab_file_RESULT.pdf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_count_mutations_script, "--gtf", str_count_mutations_gtf_file, "--key", str_count_mutations_key_genes,
                                  "--tab", str_count_mutations_input_tab_file, "--out_file", str_count_mutations_output_file_pdf ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_count_mutations_output_file_answer, str_count_mutations_output_file_result )
        # Destroy environment
        self.func_remove_files( [ str_count_mutations_output_file_result, str_count_mutations_output_file_pdf ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def not_test_tabs_to_percent_mutations_for_one_tab_file_secondary( self ):
        """
        Test tabs_to_percent_mutations.py for reading in one tab file with mutations, using secondary evidence.
        """
        # Create test environment
        str_count_mutations_script = os.path.join( self.str_script_dir, "tabs_to_percent_mutations.py" )
        str_count_mutations_gtf_file = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations.gtf" )
        str_count_mutations_key_genes = "Gene1,Gene3,Gene5"
        str_count_mutations_input_tab_file = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_1.tab" )
        str_count_mutations_output_file_answer = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_for_one_tab_file_secondary_ANSWER.txt" )
        str_count_mutations_output_file_result = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_one_tab_file_secondary_RESULT.txt" )
        str_count_mutations_output_file_pdf = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_one_tab_file_secondary_RESULT.pdf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_count_mutations_script, "--gtf", str_count_mutations_gtf_file, "--key", str_count_mutations_key_genes,
                                  "--tab", str_count_mutations_input_tab_file, "--out_file", str_count_mutations_output_file_pdf, "--second" ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_count_mutations_output_file_answer, str_count_mutations_output_file_result )
        # Destroy environment
        self.func_remove_files( [ str_count_mutations_output_file_result, str_count_mutations_output_file_pdf ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def not_test_tabs_to_percent_mutations_for_three_tab_file( self ):
        """
        Test tabs_to_percent_mutations.py for reading in three tab files with mutations.
        """
        # Create test environment
        str_count_mutations_script = os.path.join( self.str_script_dir, "tabs_to_percent_mutations.py" )
        str_count_mutations_gtf_file = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations.gtf" )
        str_count_mutations_key_genes = "Gene1,Gene3,Gene5"
        str_count_mutations_input_tab_file_1 = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_1.tab" )
        str_count_mutations_input_tab_file_2 = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_2.tab" )
        str_count_mutations_input_tab_file_3 = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_3.tab" )
        str_count_mutations_output_file_answer = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_for_three_tab_file_ANSWER.txt" )
        str_count_mutations_output_file_result = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_three_tab_file_RESULT.txt" )
        str_count_mutations_output_file_pdf = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_three_tab_file_RESULT.pdf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_count_mutations_script, "--gtf", str_count_mutations_gtf_file, "--key", str_count_mutations_key_genes,
                                  "--tab", str_count_mutations_input_tab_file_1, "--tab", str_count_mutations_input_tab_file_2, "--tab", str_count_mutations_input_tab_file_3,
                                  "--out_file", str_count_mutations_output_file_pdf ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_count_mutations_output_file_answer, str_count_mutations_output_file_result )
        # Destroy environment
        self.func_remove_files( [ str_count_mutations_output_file_result, str_count_mutations_output_file_pdf ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def not_test_tabs_to_percent_mutations_for_three_tab_file_secondary( self ):
        """
        Test tabs_to_percent_mutations.py for reading in three tab files with mutations, using secondary evidence.
        """
        # Create test environment
        str_count_mutations_script = os.path.join( self.str_script_dir, "tabs_to_percent_mutations.py" )
        str_count_mutations_gtf_file = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations.gtf" )
        str_count_mutations_key_genes = "Gene1,Gene3,Gene5"
        str_count_mutations_input_tab_file_1 = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_1.tab" )
        str_count_mutations_input_tab_file_2 = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_2.tab" )
        str_count_mutations_input_tab_file_3 = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_tab_3.tab" )
        str_count_mutations_output_file_answer = os.path.join( self.str_test_data_dir, "tabs_to_percent_mutations_for_three_tab_file_secondary_ANSWER.txt" )
        str_count_mutations_output_file_result = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_three_tab_file_secondary_RESULT.txt" )
        str_count_mutations_output_file_pdf = os.path.join( self.str_test_data_dir_working, "tabs_to_percent_mutations_for_three_tab_file_secondary_RESULT.pdf" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_count_mutations_script, "--gtf", str_count_mutations_gtf_file, "--key", str_count_mutations_key_genes,
                                  "--tab", str_count_mutations_input_tab_file_1, "--tab", str_count_mutations_input_tab_file_2, "--tab", str_count_mutations_input_tab_file_3,
                                  "--out_file", str_count_mutations_output_file_pdf, "--second" ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_count_mutations_output_file_answer, str_count_mutations_output_file_result )
        # Destroy environment
        self.func_remove_files( [ str_count_mutations_output_file_result, str_count_mutations_output_file_pdf ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def not_test_vcfs_to_snp_calls_tab_filter_maf_vcf( self ):
        """
        Test vcfs_to_snp_calls_tab.py with filtering. Inputs are maf and vcf files.
        """
        # Create test environment
        str_snp_calls_script = os.path.join( self.str_script_dir, "vcfs_to_snp_calls_tab.py" )
        str_snp_calls_input_file_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_filter.maf" )
        str_maf_tumor_key = "test"
        str_snp_calls_input_file_2 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.vcf" )
        str_snp_calls_input_depth_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.depth" )
        str_snp_calls_input_depth_2 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.depth" )
        str_snp_calls_answer = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_filter_maf_vcf_1_ANSWER_sorted.tab" )
        str_snp_calls_result = os.path.join( self.str_test_data_dir_working, "vcfs_to_snp_calls_tab_filter_maf_vcf_1_RESULT.tab" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_snp_calls_script, "--maf_reference", str_snp_calls_input_file_1, "--tumor", str_maf_tumor_key,
                                  "--vcf", str_snp_calls_input_file_2, "--count_reference", str_snp_calls_input_depth_1,
                                  "--count", str_snp_calls_input_depth_2, str_snp_calls_result ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        lstr_answer_lines = None
        with open( str_snp_calls_answer, "r" ) as hndl_answer:
            lstr_answer_lines = [ str_line for str_line in hndl_answer.read().split("\n") if str_line ]
        lstr_answer_lines.sort()
        lstr_result_lines = None
        with open( str_snp_calls_result, "r" ) as hndl_result:
            lstr_result_lines = [ str_line for str_line in hndl_result.read().split("\n") if str_line ]
        lstr_result_lines.sort()
        # Destroy environment
        self.func_remove_files( [ str_snp_calls_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_equals( "\n".join( lstr_answer_lines), "\n".join( lstr_result_lines ) )

    def not_test_vcfs_to_snp_calls_tab_filter_maf_vcf_2( self ):
        """
        Test vcfs_to_snp_calls_tab.py with filtering. Inputs are maf and vcf 2 files.
        """
        # Create test environment
        str_snp_calls_script = os.path.join( self.str_script_dir, "vcfs_to_snp_calls_tab.py" )
        str_snp_calls_input_file_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_filter.maf" )
        str_maf_tumor_key = "test"
        str_snp_calls_input_file_2 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_2_filter.vcf" )
        str_snp_calls_input_depth_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.depth" )
        str_snp_calls_input_depth_2 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_2_filter.depth" )
        str_snp_calls_answer = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_filter_maf_vcf_2_ANSWER_sorted.tab" )
        str_snp_calls_result = os.path.join( self.str_test_data_dir_working, "vcfs_to_snp_calls_tab_filter_maf_vcf_2_RESULT.tab" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_snp_calls_script, "--maf_reference", str_snp_calls_input_file_1, "--tumor", str_maf_tumor_key,
                                  "--vcf", str_snp_calls_input_file_2, "--count_reference", str_snp_calls_input_depth_1,
                                  "--count", str_snp_calls_input_depth_2, str_snp_calls_result ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        lstr_answer_lines = None
        with open( str_snp_calls_answer, "r" ) as hndl_answer:
            lstr_answer_lines = [ str_line for str_line in hndl_answer.read().split("\n") if str_line ]
        lstr_answer_lines.sort()
        lstr_result_lines = None
        with open( str_snp_calls_result, "r" ) as hndl_result:
            lstr_result_lines = [ str_line for str_line in hndl_result.read().split("\n") if str_line ]
        lstr_result_lines.sort()
        # Destroy environment
        self.func_remove_files( [ str_snp_calls_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_equals( "\n".join( lstr_answer_lines), "\n".join( lstr_result_lines ) )

    def not_test_vcfs_to_snp_calls_tab_filter_vcf_1_2( self ):
        """
        Test vcfs_to_snp_calls_tab.py with filtering. Inputs are and vcf 1 and 2 files.
        """
        # Create test environment
        str_snp_calls_script = os.path.join( self.str_script_dir, "vcfs_to_snp_calls_tab.py" )
        str_snp_calls_input_file_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.vcf" )
        str_snp_calls_input_file_2 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_2_filter.vcf" )
        str_snp_calls_input_depth_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.depth" )
        str_snp_calls_input_depth_2 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_2_filter.depth" )
        str_snp_calls_answer = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_filter_vcf_1_2_ANSWER_sorted.tab" )
        str_snp_calls_result = os.path.join( self.str_test_data_dir_working, "vcfs_to_snp_calls_tab_filter_vcf_1_2_RESULT.tab" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_snp_calls_script, "--vcf_reference", str_snp_calls_input_file_1,
                                  "--vcf", str_snp_calls_input_file_2, "--count_reference", str_snp_calls_input_depth_1,
                                  "--count", str_snp_calls_input_depth_2, str_snp_calls_result ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        lstr_answer_lines = None
        with open( str_snp_calls_answer, "r" ) as hndl_answer:
            lstr_answer_lines = [ str_line for str_line in hndl_answer.read().split("\n") if str_line ]
        lstr_answer_lines.sort()
        lstr_result_lines = None
        with open( str_snp_calls_result, "r" ) as hndl_result:
            lstr_result_lines = [ str_line for str_line in hndl_result.read().split("\n") if str_line ]
        lstr_result_lines.sort()
        # Destroy environment
        self.func_remove_files( [ str_snp_calls_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_equals( "\n".join( lstr_answer_lines), "\n".join( lstr_result_lines ) )

    def not_test_vcfs_to_genotype_matrix_1_file( self ):
        """
        Test vcfs_to_genotype_matrix.py with one input file.
        """
        # Create test environment
        str_genotype_script = os.path.join( self.str_script_dir, "vcfs_to_genotype_matrix.py" )
        str_vcf_directory = os.path.join( self.str_test_data_dir, "test_vcf_genotype_matrix_dir_1" )
        str_genotype_answer = os.path.join( self.str_test_data_dir, "vcfs_to_genotype_matrix_1_ANSWER.txt" )
        str_genotype_result = os.path.join( self.str_test_data_dir_working, "vcfs_to_genotype_matrix_1_RESULT.txt" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_genotype_script, "--matrix", str_genotype_result, str_vcf_directory ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_genotype_answer, str_genotype_result )
        # Destroy environment
        self.func_remove_files( [ str_genotype_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def not_test_vcfs_to_genotype_matrix_3_file( self ):
        """
        Test vcfs_to_genotype_matrix.py with one input file in one directory and 2 in another.
        """
        # Create test environment
        str_genotype_script = os.path.join( self.str_script_dir, "vcfs_to_genotype_matrix.py" )
        str_vcf_directory_1 = os.path.join( self.str_test_data_dir, "test_vcf_genotype_matrix_dir_1" )
        str_vcf_directory_2 = os.path.join( self.str_test_data_dir, "test_vcf_genotype_matrix_dir_2" )
        str_genotype_answer = os.path.join( self.str_test_data_dir, "vcfs_to_genotype_matrix_3_ANSWER.txt" )
        str_genotype_result = os.path.join( self.str_test_data_dir_working, "vcfs_to_genotype_matrix_3_RESULT.txt" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call Example script
        str_command = " ".join( [ str_genotype_script, "--matrix", str_genotype_result, str_vcf_directory_1, str_vcf_directory_2 ] )
        Commandline.Commandline().func_CMD( str_command )
        # Check test environment for results
        f_success = self.func_are_files_equivalent( str_genotype_answer, str_genotype_result )
        # Destroy environment
        self.func_remove_files( [ str_genotype_result ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def not_test_visualize_mutation_depth_tab_files_for_error_counts_opt( self ):
        """
        Tests to make sure the TP, FP, FN, senstivity, and specificity measurements are correct from a test data set.
        This is testing output that has a changing feature space (optimization figure) and not the "ROC" plot.
        """
        # Create environment
        str_vis_script = os.path.join( self.str_script_dir, "visualize_mutation_depth_tab_files.R" )
        str_test_input_file = os.path.join( self.str_test_data_dir, "test_visualize_tab.tab" )
        str_answer_file = os.path.join( self.str_test_data_dir, "test_visualize_mutation_depth_tab_files_for_error_counts_opt_ANSWER.txt" )
        str_result_file = os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_data.txt" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call example script
        str_command = " ".join( [ str_vis_script, "-o", self.str_test_data_dir_working, "-k RNA_DNA", str_test_input_file ])
        Commandline.Commandline().func_CMD( str_command )
        # Check for sucess
        f_success = self.func_are_files_equivalent( str_answer_file, str_result_file )
        # Destroy environment
        self.func_remove_files( [ str_result_file, os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_depth_distributions.pdf" ),
                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_fdr_min_read_coverage_norm.pdf" ),
                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_data_truth_held.txt" ),
                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_optimize_detail_validation.pdf" ),
                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_raw_class_distributions_detail_validation.pdf" ),
                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_roc_detail_validation.pdf" ),
                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab.tab_sensitivity_min_read_coverage_norm.pdf" ) ] )
        self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success )

    def test_visualize_mutation_depth_tab_files_for_error_counts_roc( self ):
        """
        Tests to make sure the TP, FP, FN, senstivity, and specificity measurements are correct from a test data set.
        This is testing output that has a set feature space (the "ROC" plot) and not the optimization plot.
        """
        # Create environment
        str_vis_script = os.path.join( self.str_script_dir, "visualize_mutation_depth_tab_files.R" )
        str_test_input_file = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab" )
        str_answer_file_1 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_1_answer.txt" )
        str_answer_file_2 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_2_answer.txt" )
        str_answer_file_3 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_3_answer.txt" )
        str_answer_file_4 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_4_answer.txt" )
        str_result_file_1 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_1.txt" )
        str_result_file_2 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_2.txt" )
        str_result_file_3 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_3.txt" )
        str_result_file_4 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_4.txt" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call example script
        str_command = " ".join( [ str_vis_script, "-o", self.str_test_data_dir_working, "-k RNA_DNA", str_test_input_file ])
        Commandline.Commandline().func_CMD( str_command )
        # Check for sucess
        f_success_1 = self.func_are_files_equivalent( str_answer_file_1, str_result_file_1 )
        f_success_2 = self.func_are_files_equivalent( str_answer_file_2, str_result_file_2 )
        f_success_3 = self.func_are_files_equivalent( str_answer_file_3, str_result_file_3 )
        f_success_4 = self.func_are_files_equivalent( str_answer_file_4, str_result_file_4 )
        # Destroy environment
#        self.func_remove_files( [ os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_depth_distributions.pdf" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_fdr_min_read_coverage_norm.pdf" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data.txt" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_raw_class_distributions_detail_validation.pdf" ),
 #                                 os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_roc.pdf" ),
 #                                 os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_optimize_detail_validation.pdf" ),
 #                                 os.path.join( str_result_file_1 ),
 #                                 os.path.join( str_result_file_2 ),
 #                                 os.path.join( str_result_file_3 ),
 #                                 os.path.join( str_result_file_4 ),
 #                                 os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_sensitivity_min_read_coverage_norm.pdf" ) ] )
 #       self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success_1 and f_success_2 and f_success_3 and f_success_4 )

    def not_test_visualize_mutation_depth_tab_files_for_roc_like_rnaseq( self ):
        """
        Tests to make sure the TP, FP, FN, senstivity, and specificity measurements are correct from a test data set.
        This is testing output that has a set feature space (the "ROC" plot) and not the optimization plot.

        The other test uses a simple input data set similar to traditional ROC data, this one have varying RNA seq depth and 
        such that allows a more authentic test.
        """
        # Create environment
        str_vis_script = os.path.join( self.str_script_dir, "visualize_mutation_depth_tab_files.R" )
        str_test_input_file = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc_like_rnaseq.tab" )
        str_answer_file_1 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_1_answer.txt" )
        str_answer_file_2 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_like_ranseq_2_answer.txt" )
        str_answer_file_3 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_3_answer.txt" )
        str_answer_file_4 = os.path.join( self.str_test_data_dir, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_4_answer.txt" )
        str_result_file_1 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_1.txt" )
        str_result_file_2 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_2.txt" )
        str_result_file_3 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_3.txt" )
        str_result_file_4 = os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc.tab_data_roc_like_rnaseq_4.txt" )
        self.func_make_dummy_dir( self.str_test_data_dir_working )
        # Call example script
        str_command = " ".join( [ str_vis_script, "-o", self.str_test_data_dir_working, "-k RNA_DNA", str_test_input_file ])
        Commandline.Commandline().func_CMD( str_command )
        # Check for sucess
        f_success_1 = self.func_are_files_equivalent( str_answer_file_1, str_result_file_1 )
        f_success_2 = self.func_are_files_equivalent( str_answer_file_2, str_result_file_2 )
        f_success_3 = self.func_are_files_equivalent( str_answer_file_3, str_result_file_3 )
        f_success_4 = self.func_are_files_equivalent( str_answer_file_4, str_result_file_4 )
        # Destroy environment
#        self.func_remove_files( [ os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_depth_distributions.pdf" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_fdr_min_read_coverage_norm.pdf" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_data.txt" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_raw_class_distributions_detail_validation.pdf" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_roc.pdf" ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_optimize_detail_validation.pdf" ),
#                                  os.path.join( str_result_file_1 ),
#                                  os.path.join( str_result_file_2 ),
#                                  os.path.join( str_result_file_3 ),
#                                  os.path.join( str_result_file_4 ),
#                                  os.path.join( self.str_test_data_dir_working, "test_visualize_tab_roc_like_rnaseq.tab_sensitivity_min_read_coverage_norm.pdf" ) ] )
 #       self.func_remove_dirs( [ self.str_test_data_dir_working ] )
        # Evaluate
        self.func_test_true( f_success_1 and f_success_2 and f_success_3 and f_success_4 )

# Creates a suite of tests
def suite():
    return unittest.TestLoader().loadTestsFromTestCase( ScriptTester )
