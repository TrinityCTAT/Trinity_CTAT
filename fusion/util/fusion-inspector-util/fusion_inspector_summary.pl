#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 samples.txt processing_dir\n\n";

my $samples_file = $ARGV[0] or die $usage;
my $processing_dir = $ARGV[1] or die $usage;

main: {

    my %samples;
    {
        open (my $fh, $samples_file) or die "Error, cannot open file $samples_file";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $sample = $x[0];
            $samples{$sample} = 1;
        }
        close $fh;
        
    }

    foreach my $sample (sort keys %samples) {
        my $summary_report = "$processing_dir/samples/$sample/CHIM_VISUALIZER/chim_vis.fusion_preds.coalesced.summary.abridged";
        open (my $fh, $summary_report) or die "Error, cannot open file $summary_report";
        while (<$fh>) {
            my @x = split(/\t/);
            my ($geneA, $geneB) = ($x[0], $x[2]);
            
            print "$sample|$geneA--$geneB\t" . join("\t", @x);
        }
        close $fh;
    }
    
    exit(0);
}

