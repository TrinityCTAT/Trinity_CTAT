#!/usr/bin/env python

import sys, os, re
import subprocess

usage = "usage: {} input.bam output.bam\n\n".format(sys.argv[0])

if len(sys.argv) < 3:
    sys.stderr.write(usage)
    sys.exit(1)

input_bam_file = sys.argv[1]
output_bam_file = sys.argv[2]

    
def main():

    if reads_need_cleaning(input_bam_file):
        sys.stderr.write("cleaning reads\n")
        clean_reads(input_bam_file, output_bam_file)
    else:
        sys.stderr.write("reads look fine, just symlinking\n")
        subprocess.check_call("ln -s {} {}".format(input_bam_file, output_bam_file), shell=True)

    sys.exit(0)


def reads_need_cleaning(bam_file):

    bam_reader = subprocess.Popen(['samtools', 'view', bam_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    counter = 0
    for line in bam_reader.stdout:
        fields = line.split("\t")
        readname = fields[0]
        counter += 1
        if re.search("/[12]$", readname):
            return True

        if counter >= 100:
            break
        
    return False

def clean_reads(input_bam, output_bam):

    sam_reader = subprocess.Popen(["samtools", "view", "-h", input_bam], stdout=subprocess.PIPE)
    sam_writer = subprocess.Popen(["samtools", "view", "-Sb", "-o", output_bam], stdin=subprocess.PIPE)

    for line in sam_reader.stdout:
        x = line.split("\t")
        x[0] = re.sub("/[12]$", "", x[0])
        sam_writer.stdin.write("\t".join(x))
    


if __name__ == '__main__':
    main()
