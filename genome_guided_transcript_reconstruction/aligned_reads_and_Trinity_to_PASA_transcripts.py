#!/usr/bin/env python


__author__ = "Brian Haas"
__copyright__ = "Copyright 2015"
__credits__ = [ "Timothy Tickle", "Brian Haas" ]
__license__ = "BSD 3-clause"
__maintainer__ = "Brian Haas"
__email__ = "bhaas@broadinstitute.org"
__status__ = "Development"

import inspect
import os,sys
import sciedpiper.Command as Command
import sciedpiper.ParentScript as ParentScript


FUSION_ANNOTATOR_LIB=""
if not os.environ.has_key("FUSION_ANNOTATOR_LIB"):
    print >> sys.stderr, "Error, need FUSION_ANNOTATOR_LIB environmental variable set and pointing to FUSION_ANNOTATOR_LIB installation directory"
    sys.exit(2)
else:
    FUSION_ANNOTATOR_LIB = os.environ["FUSION_ANNOTATOR_LIB"]


PASA_LITE_HOME=""
if not os.environ.has_key("PASA_LITE_HOME"):
    print >> sys.stderr, "Error, need PASA_LITE_HOME environmental variable set to the PASA_LITE installation directory"
    sys.exit(2)
else:
    PASA_LITE_HOME = os.environ["PASA_LITE_HOME"]


class TranscriptReconstructionPipeline( ParentScript.ParentScript ):
    
    def func_update_arguments(self, arg_raw ):
        """
        Updates to the arg parser, command line options
        
        * arg_raw : Arguments ( not yet parsed )
                  : Arguments
        * return  : Updated Arguments
                  : Arguments
        """

        arg_raw.prog = "aligned_reads_and_Trinity_to_PASA_transcripts.py"
        arg_raw.description = str("assembles transcripts from read alignments using cufflinks, " +
                                  "then uses PASA_Lite to assemble cufflinks transcripts with Trinity alignments")
        arg_raw.add_argument("--bam", dest = "bam_file", required=True, default = None, help = "coordinate-sorted bam file" )        
        arg_raw.add_argument("--trinity_gff3_file", required=True, dest="trinity_gff3_file", default=None, help = "gff3 file containing gmap-aligned Trinity transcripts")
        arg_raw.add_argument("--CPU", dest="CPU", default=1, help="number of threads to use")

    def func_make_commands( self, args_parsed, cur_pipeline ):
        """
        Allows:
        - the creation of commands in the child object.
        - the creation of directories.
        - checking that files exist.
        
        To know the variables available from command line look in the ParentScript in func_create_arguments.
        """


        CPU = args_parsed.CPU
        
                
        # Make commands
        lcmd_commands = []
        
        ################
        # run cufflinks
        ################
        
        bam_file = args_parsed.bam_file
        cufflinks_gtf_file = "cufflinks_outdir/transcripts.gtf"
        
        cmdstr = str("cufflinks " +
                     # " -g FUSION_ANNOTATOR_LIB + "/gencode.v19.rna_seq_pipeline.gtf " +
                     " -o cufflinks_outdir " + bam_file)
        
        
        lcmd_commands.append( Command.Command( str_cur_command = cmdstr,
                                               lstr_cur_dependencies = [ bam_file ], 
                                               lstr_cur_products = [ cufflinks_gtf_file ] ) )
        

        ###################
        ## PASA Lite ######
        ###################

        trinity_gff3_file = args_parsed.trinity_gff3_file
        pasa_lite_valid_alignments = "pasa_lite.valid_alignments.gtf"
        hg19_fa = FUSION_ANNOTATOR_LIB + "/Hg19.fa"

        ## Alignment validation
        cmdstr = str(PASA_LITE_HOME + "/PASA.alignmentValidator --genome " + hg19_fa + " --CPU " + `CPU` + " " +
                     trinity_gff3_file + " " + cufflinks_gtf_file)

        lcmd_commands.append( Command.Command( str_cur_command = cmdstr,
                                               lstr_cur_dependencies = [ trinity_gff3_file, cufflinks_gtf_file, hg19_fa ],
                                               lstr_cur_products = [ pasa_lite_valid_alignments ] ) )


        ## Alignment assembly
        cmdstr = str(PASA_LITE_HOME + "/PASA.alignmentAssembler --CPU " + `CPU` + " " + pasa_lite_valid_alignments)

        pasa_assemblies_file = "pasa_lite.pasa_assembled_alignments.gtf"
        
        lcmd_commands.append( Command.Command( str_cur_command = cmdstr,
                                               lstr_cur_dependencies = [ pasa_lite_valid_alignments ],
                                               lstr_cur_products = [ pasa_assemblies_file ] ) )
                    

        ## Generate fasta sequences
        pasa_fasta_output_filename = "pasa_lite.pasa_assembled_alignments.fasta"
        
        cmdstr = str(PASA_LITE_HOME + "/util/gfx_alignment_to_cdna_fasta.pl " + pasa_assemblies_file + " " +
                     hg19_fa + " > " + pasa_fasta_output_filename)

        lcmd_commands.append( Command.Command( str_cur_command = cmdstr,
                                               lstr_cur_dependencies = [ pasa_assemblies_file, hg19_fa ],
                                               lstr_cur_products = [ pasa_fasta_output_filename ] ) )
        
        
        return lcmd_commands
    
    
if __name__ == "__main__":

    # Needed to run, calls the script
    TranscriptReconstructionPipeline().func_run_pipeline()
