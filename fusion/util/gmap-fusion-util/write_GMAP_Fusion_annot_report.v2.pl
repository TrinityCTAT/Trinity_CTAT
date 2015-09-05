#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");

use __GLOBALS__;
use FusionAnnotator;

my $usage = "usage: $0 gmap.map.gff3.chims_described.w_read_support\n\n";

my $chim_file = $ARGV[0] or die $usage;

main: {

    open (my $fh, $chim_file) or die "Error, cannot open file $chim_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $trinity_acc = $x[0];
        my $fusion_info = $x[3];
        my $J = $x[4];
        my $S = $x[5];
        
        if (! defined($J)) { $J = "."; }
        if (! defined($S)) { $S = "."; }

        my ($geneA, $deltaA, $trans_brkptA, $breakpointA, 
            $geneB, $deltaB, $trans_brkptB, $breakpointB, $fusion_name) = split(/;/, $fusion_info);
        
        
        my ($chrA, $coordA) = split(/:/, $breakpointA);
        my ($chrB, $coordB) = split(/:/, $breakpointB);
        
        my $trans_brkpt = join("-", $trans_brkptA, $trans_brkptB);
        
        my ($junction_type) = ($deltaA != 0 || $deltaB != 0) ? "INCL_NON_REF_SPLICE" : "ONLY_REF_SPLICE";
        
        print join("\t", $fusion_name, $trinity_acc, $trans_brkpt, 
                   $geneA, $chrA, $coordA, $geneB, $chrB, $coordB, 
                   $J, $S, $junction_type) . "\n";
    }
    
    close $fh;

    exit(0);
    
}
