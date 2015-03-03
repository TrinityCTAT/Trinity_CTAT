
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

# Creates a suite of tests
def suite():
    return unittest.TestLoader().loadTestsFromTestCase( ScriptTester )
