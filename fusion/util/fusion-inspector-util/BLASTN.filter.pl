#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);                                                 



####
# All this mainly does is to reformat the fusion inspector predictions output so that we can filter it using the STAR-Fusion.filter script.
# However, also filtering out those non-consensus splice types having less than 10 junction reads supporting them.
####


######################################################

my $STAR_FUSION_DIR = $ENV{STAR_FUSION_DIR} or die "Error, need env var for STAR_FUSION_DIR";


my $usage = <<__EOUSAGE__;

########################################################################
#
# Required:
#
#  --fusion_preds <string>        preliminary fusion predictions
#                                 Required formatting is:  
#                                 geneA--geneB (tab) score (tab) ... rest
#
#    --min_novel_junction_support <int>    default: 10  (minimum of 10 junction reads required if breakpoint
#                                                        lacks involvement of only reference junctions)
#
#    --min_alt_pct_junction <float>        default: 10.0  (10% of the dominant isoform junction support)
#
########################################################################  


__EOUSAGE__

    ;

my $help_flag;

my $fusion_preds_file;

my $min_novel_junction_support = 10;
my $min_alt_pct_junction = 10.0;


&GetOptions ( 'h' => \$help_flag, 
              
              'fusion_preds=s' => \$fusion_preds_file,
              
              'min_novel_junction_support=i' => \$min_novel_junction_support,
              'min_alt_pct_junction=f' => \$min_alt_pct_junction,
    );


if ($help_flag) {
    die $usage;
}

unless ($fusion_preds_file) {
    die $usage;
}


main: {

    my $star_fusion_fmt_file = "$fusion_preds_file.starFfmt";
    open (my $ofh, ">$star_fusion_fmt_file") or die  "Error, cannot write to $star_fusion_fmt_file";
    
    
    my @fusions;
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    my $header = <$fh>;
    
    
    
    my @fusions;

    while (<$fh>) {
        if (/^\#/) { 
            next; 
        }
        chomp;
        my $line = $_;

        my ($geneA, $chr_brkpt_A, $geneB, $chr_brkpt_B, $splice_type, $junction_count, $spanning_count, $num_left_contrary_reads, $num_right_contrary_reads, $TAF_left, $TAF_right, $fusion_annotations) = split(/\t/);

        my $fusion_name = "$geneA--$geneB";
        
        my @data = ($fusion_name, $junction_count, $spanning_count, $splice_type, $geneA, $chr_brkpt_A, $geneB, $chr_brkpt_B, $num_left_contrary_reads, $TAF_left, $num_right_contrary_reads, $TAF_right, $fusion_annotations);
        
        push (@fusions, [@data]);

        #print $ofh join("\t", 
        
    }
    

    ## filter and print
    print $ofh join("\t", "#fusion_name", "junction_count", "spanning_count", "splice_type", "geneA", "chr_brktp_A", "geneB", "chr_brktp_B", "num_left_contrary", "num_right_contrary", "TAF_left", "TAF_right", "fusion_annotations") . "\n";

    @fusions = reverse sort {$a->[1]<=>$b->[1]} @fusions;


    my %dominant_junction_support;

    foreach my $fusion (@fusions) {
        
        my $fusion_name = $fusion->[0];
        my $splice_type = $fusion->[3];
        my $junction_support = $fusion->[1];
        if ($splice_type eq 'INCL_NON_REF_SPLICE' && $junction_support < $min_novel_junction_support) {
            next;
        }
        if (! exists $dominant_junction_support{$fusion_name}) {
            $dominant_junction_support{$fusion_name} = $junction_support;
        }
        else {
            if ($junction_support < $dominant_junction_support{$fusion_name} * $min_alt_pct_junction/100) {
                next;
            }
        }

        print $ofh join("\t", @$fusion) . "\n";
    }
            
    close $fh;
    close $ofh;

    ## now do the homology filter
    
    my $cmd = "$STAR_FUSION_DIR/util/STAR-Fusion.filter --fusion_preds $star_fusion_fmt_file";
    &process_cmd($cmd);

}



####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
        
    my $ret = system($cmd);
    if ($ret) {

        die "Error, cmd $cmd died with ret $ret";
    }
    
    return;
}
    
        
