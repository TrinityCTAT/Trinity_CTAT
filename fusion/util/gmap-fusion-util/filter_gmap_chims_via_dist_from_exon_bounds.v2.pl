#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 gmap.gff3.summary minJ minJS min_novel_J\n\n";

my $gmap_summary_file = $ARGV[0] or die $usage;

my $minJ = $ARGV[1];
my $minJS = $ARGV[2];
my $min_novel_J = $ARGV[3];

unless (defined ($minJ) && defined($minJS) && defined($min_novel_J)) {
    die $usage;
}


open (my $fh, $gmap_summary_file) or die $!;
while (<$fh>) {
    if (/^\#/) { # header
        print $_;
        next;
    }

    chomp;
    my $line = $_;
    my @x = split(/\t/);
    my $info = $x[3];
    
    my $J = $x[4];
    my $S = $x[5];
    
    
    if ($info ne '.') {
        my ($geneA, $trans_brktpA, $distA, $breakpointA, 
            $geneB, $trans_brkptB, $distB, $breakpointB, 
            $chim) = split(/;/, $info);

        if ($geneA eq $geneB) { next; }
        
        my $record_pass = 0;

        if (defined($J) && defined($S) && $J =~ /^\d+$/ && $S =~ /^\d+$/) {
            
            ## Examining Illumina PE read support

            my $sumJS = $J + $S;
            if ($distA == 0 && $distB == 0 && $J >= $minJ && $sumJS >= $minJS) {
                # reference junction
                $record_pass = 1;
            }
            elsif ( ($distA > 0 || $distB > 0) && $J >= $min_novel_J  && $sumJS >= $minJS) {
                $record_pass = 1;
            }
        
        }
        else {
            ## only reporting ref-junction entries
            
            if ($distA == 0 && $distB == 0) {
                $record_pass = 1;
            }
        }
        
        if ($record_pass) {
            print "$line\t$chim\n";
        }
    }
    
}


exit(0);
