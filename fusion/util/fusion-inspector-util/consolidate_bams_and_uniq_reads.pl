#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\tusage: $0 target.fasta fileA.bam,fileB.bam,...  out_prefix\n\n";

my $target_fa = $ARGV[0] or die $usage;
my $bam_files = $ARGV[1] or die $usage;
my $out_prefix = $ARGV[2] or die $usage;


main: {

    my @bams = split(/,/, $bam_files);

    my @tmp_files;

    my $tmp_sam = "$out_prefix.tmp.sam";
    if (-s $tmp_sam) {
        unlink($tmp_sam) or die "Error, cannot remove tmp file $tmp_sam";
    }
    
    foreach my $bam (@bams) {
        
        my $cmd = "samtools view $bam >> $tmp_sam";
        &process_cmd($cmd);
    }
    push (@tmp_files, $tmp_sam);
    

    # sort by scaff, readname, coordinate
    my $cmd = "sort -k3,3 -k1,1 -k4,4n $tmp_sam > $tmp_sam.nameSorted";
    &process_cmd($cmd);
    push (@tmp_files, "$tmp_sam.nameSorted");
    
    ## uniq it
    &uniq_the_reads("$tmp_sam.nameSorted", "$tmp_sam.nameSorted.uniq");
    push (@tmp_files, "$tmp_sam.nameSorted.uniq");

    ## convert to coordinate sorted bam
    $cmd = "bash -c \"set -o pipefail && cat $tmp_sam.nameSorted.uniq | samtools view -Sb -T $target_fa - | samtools sort - $out_prefix.cSorted \"";
    &process_cmd($cmd);
    
    $cmd = "samtools index $out_prefix.cSorted.bam";
    &process_cmd($cmd);

    unlink(@tmp_files);
    

    exit(0);
}

####
sub uniq_the_reads {
    my ($in_sam, $out_sam) = @_;

    my %seen;
    my $prev_core_read_name = "";

    open (my $ofh, ">$out_sam") or die "Error, cannot write to $out_sam";

    my $sam_reader = new SAM_reader($in_sam);
    while (my $sam_entry = $sam_reader->get_next()) {



        my $full_read_name = $sam_entry->reconstruct_full_read_name();
        my $core_read_name = $sam_entry->get_core_read_name();

        if ($core_read_name ne $prev_core_read_name) {
            %seen = ();
            $prev_core_read_name = $core_read_name;
        }

        unless ($seen{$full_read_name}) {
            $seen{$full_read_name} = 1;
            print $ofh $sam_entry->get_original_line() . "\n";
        }
    }
    close $ofh;

    return;
}
    
####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}
        
