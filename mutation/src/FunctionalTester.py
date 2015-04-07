
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

class FunctionalTester( ParentPipelineTester.ParentPipelineTester ):
    """
    Functional testing for scripts, these are focused on making sure the script can be called in different ways and complete without error.
    """
    
    # Testing environment
    str_script_dir = "/ahg/regev/users/ttickle/dev/KCO/SOFTWARE/Trinity_CTAT/mutation/src"
    str_test_data = "/ahg/regev/users/ttickle/dev/KCO/SOFTWARE/Trinity_CTAT/mutation/demo_data"
    str_testing_area = "/ahg/regev/users/ttickle/dev/KCO/SOFTWARE/Trinity_CTAT/mutation/active_testing_script_tester"

    str_input_index = os.path.join( str_test_data, "star_index_test" )
    str_input_test_bam = os.path.join( str_test_data, "star_aligned_test.bam" )
    str_left_file = os.path.join( str_test_data, "Left_10000.fq" )
    str_right_file = os.path.join( str_test_data, "Right_10000.fq" )
    str_reference_vcf = os.path.join( str_test_data, "dbsnp_138.b37_Hg19_head_1000.vcf" )
    str_reference_genome = os.path.join( str_test_data, "Hg19_21_demo.fasta" )
    str_update_command = "--update AddOrReplaceReadGroups.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,MarkDuplicates.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,SortSam.jar:/seq/regev_genome_portal/SOFTWARE/Picard/current,GenomeAnalysisTK.jar:/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.1-1-g07a4bf8"

    def test_rnaseq_mutation_pipeline_for_no_args( self ):
        """
        Tests rnaseq_mutation_pipeline.py for no args call.
        """
        # Create test environment
        str_command = "python rnaseq_mutation_pipeline.py"
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertFalse( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_args_short( self ):
        """
        Tests rnaseq_mutation_pipeline.py for help args call short.
        """
        # Create test environment
        str_command = "python rnaseq_mutation_pipeline.py -h"
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_args_long( self ):
        """
        Tests rnaseq_mutation_pipeline.py for help args call short.
        """
        # Create test environment
        str_command = "python rnaseq_mutation_pipeline.py --help"
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_samtools_call( self ):
        """
        Tests rnaseq_mutation_pipeline.py for samtools call.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_samtools"]),"--vcf",self.str_reference_vcf, self.str_update_command ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_test( self ):
        """
        Tests rnaseq_mutation_pipeline.py for test mode.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,"--test",
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_test"]),"--vcf",self.str_reference_vcf, self.str_update_command ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_gatk_call( self ):
        """
        Tests rnaseq_mutation_pipeline.py for gatk call.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode GATK --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_gatk"]),"--vcf",self.str_reference_vcf, self.str_update_command ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_compression( self ):
        """
        Tests rnaseq_mutation_pipeline.py for compression.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_compression"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--compress","archive"  ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_clean( self ):
        """
        Tests rnaseq_mutation_pipeline.py for clean.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_clean"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--clean" ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_archive( self ):
        """
        Tests rnaseq_mutation_pipeline.py for archive.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_archive"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--copy", "hold_test_runs" ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_realign( self ):
        """
        Tests rnaseq_mutation_pipeline.py for realignment.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_realign"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--realign" ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_starting_with_bam( self ):
        """
        Tests rnaseq_mutation_pipeline.py for starting with a bam.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_bam"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--bam", self.str_input_test_bam ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_star_limited( self ):
        """
        Tests rnaseq_mutation_pipeline.py for starting with start limited mode
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode LIMITED --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_limited"]),"--vcf",self.str_reference_vcf, self.str_update_command
                                ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_named_log_file( self ):
        """
        Tests rnaseq_mutation_pipeline.py for starting with a named log file.
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_samtools"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--log", os.path.join( "_".join([self.str_testing_area,"vanilla_log"]),"run.log" ) ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_starting_with_premade_index( self ):
        """
        Tests rnaseq_mutation_pipeline.py for starting with a premade index
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_premade_index"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--index", self.str_input_index ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_no_recalibration( self ):
        """
        Tests rnaseq_mutation_pipeline.py for no recalibration
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_no_recal"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--recalibrate_sam" ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_no_move( self ):
        """
        Tests rnaseq_mutation_pipeline.py for moving files
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode SAM --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_no_move"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--move", "hold_test_runs" ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )

    def test_rnaseq_mutation_pipeline_for_dnaseq( self ):
        """
        Tests rnaseq_mutation_pipeline.py for dnaseq mode
        """
        # Create test environment
        str_command = " ".join( [ "python rnaseq_mutation_pipeline.py --alignment_mode STAR --variant_call_mode GATK --threads 8 --plot --reference",
                                   self.str_reference_genome,"--left",self.str_left_file,"--right",self.str_right_file,
                                   "--out_dir","_".join([self.str_testing_area,"vanilla_dnaseq"]),"--vcf",self.str_reference_vcf, self.str_update_command,
                                   "--validate_dnaseq" ] )
        # Run command
        f_success = Commandline.Commandline().func_CMD( str_command )
        # Test error
        self.assertTrue( f_success, str_command )


# Creates a suite of tests
def suite():
    return unittest.TestLoader().loadTestsFromTestCase( FunctionalTester )
