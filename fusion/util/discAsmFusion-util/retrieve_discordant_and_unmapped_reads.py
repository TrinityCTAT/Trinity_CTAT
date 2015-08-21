#!/usr/bin/env python

import os, sys, re
import pysam
import threading

usage = "\n\n\tusage: " + sys.argv[0] + " alignments.bam left.fq right.fq\n\n"

if len(sys.argv) < 4:
    print >>sys.stderr, usage
    sys.exit(1)


alignments_bam = sys.argv[1]
left_fq_filename = sys.argv[2]
right_fq_filename = sys.argv[3]



def main():

    bam_reader = pysam.Samfile(alignments_bam)

    proper_pairs = set()

    for read in bam_reader:
        if read.is_proper_pair:
            proper_pairs.add(read.query_name)


    left_fq_extraction_thread = Unkept_read_fq_extractor(proper_pairs, left_fq_filename, left_fq_filename + ".extracted.fq")
    right_fq_extraction_thread = Unkept_read_fq_extractor(proper_pairs, right_fq_filename, right_fq_filename + ".extracted.fq")

    left_fq_extraction_thread.start()
    right_fq_extraction_thread.start()

    num_failed = 0
    for t in (left_fq_extraction_thread, right_fq_extraction_thread):
        t.join()
        if not t.success:
            num_failed += 1
            print >>sys.stderr, "Error extracting reads from file: " + t.input_fq_filename

    sys.exit(num_failed)
    

    

class Unkept_read_fq_extractor (threading.Thread):

    thread_counter = 0

    def __init__(self, keep_set, input_fq_filename, output_fq_filename):
        threading.Thread.__init__(self)

        self.keep_set = keep_set
        self.input_fq_filename = input_fq_filename
        self.output_fq_filename = output_fq_filename

        Unkept_read_fq_extractor.thread_counter += 1
        self.id = Unkept_read_fq_extractor.thread_counter

        self.success = False
        
        

    def run(self):
        ## do the extraction

        


        self.success = True


        



main()


