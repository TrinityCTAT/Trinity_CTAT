
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

    def test_reduce_vcf_to_snp_for_small_file_no_filter( self ):
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

    def test_reduce_vcf_to_snp_for_large_file_no_filter( self ):
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

    def test_reduce_vcf_to_snp_for_small_file_filter_reference( self ):
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

    def test_reduce_vcf_to_snp_for_small_file_filter( self ):
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

    def test_reduce_vcf_to_snp_for_large_file_filter( self ):
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

    def test_tabs_to_percent_mutations_for_one_tab_file( self ):
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

    def test_tabs_to_percent_mutations_for_one_tab_file_secondary( self ):
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

    def test_tabs_to_percent_mutations_for_three_tab_file( self ):
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

    def test_tabs_to_percent_mutations_for_three_tab_file_secondary( self ):
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

    def test_vcfs_to_snp_calls_tab_filter_maf_vcf( self ):
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

    def test_vcfs_to_snp_calls_tab_filter_maf_vcf_2( self ):
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

    def test_vcfs_to_snp_calls_tab_filter_vcf_1_2( self ):
        """
        Test vcfs_to_snp_calls_tab.py with filtering. Inputs are and vcf 1 and 2 files.
        """
        # Create test environment
        str_snp_calls_script = os.path.join( self.str_script_dir, "vcfs_to_snp_calls_tab.py" )
        str_snp_calls_input_file_1 = os.path.join( self.str_test_data_dir, "vcfs_to_snp_calls_tab_1_filter.vcf" )
        str_maf_tumor_key = "test"
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

# Creates a suite of tests
def suite():
    return unittest.TestLoader().loadTestsFromTestCase( ScriptTester )
