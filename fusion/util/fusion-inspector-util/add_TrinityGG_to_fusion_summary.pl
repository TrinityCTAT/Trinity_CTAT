#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "\n\tusage: $0 fusion_summary_file TrinGG_fusion_gff3\n\n";

my $fusion_summary_file = $ARGV[0] or die $usage;
my $trin_fusions_gff3 = $ARGV[1] or die $usage;


main: {

    my %breakpoint_info;
    {
        open (my $fh, $trin_fusions_gff3) or die "Error, cannot open file $trin_fusions_gff3";
        while (<$fh>) {
            chomp;
            if (/^\#TrinityFusionTranscript/) {
                chomp;
                my @x = split(/\t/);
                my ($token, $trin_acc, $contig_brkpt) = @x;

                push (@{$breakpoint_info{$contig_brkpt}}, $trin_acc);
            }
        }
        close $fh;
    }
    
    open (my $fh, $fusion_summary_file) or die "Error, cannot open file $fusion_summary_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) {
            print "$_\tTrinGG_Fusion\n";
            next;
        }
        my @x = split(/\t/);
        my ($geneA, $local_brkpt_A, $geneB, $local_brkpt_B) = ($x[0], $x[1], $x[3], $x[4]);
        my $fusion_name = "$geneA--$geneB";
        my $fusion_brkpt = "$fusion_name:$local_brkpt_A-$local_brkpt_B";
        if (my $trin_list_aref = $breakpoint_info{$fusion_brkpt}) {
            push (@x, join(",", @$trin_list_aref) );
        }
        else {
            push (@x, "."); # placeholder
        }

        print join("\t", @x) . "\n";
    }
    close $fh;

    

    exit(0);
}


            
