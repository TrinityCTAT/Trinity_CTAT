#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use __GLOBALS__;
use FusionAnnotator;

my $usage = "\n\tusage: $0 fi_test.fusion_preds.coalesced.summary.wTrinityGG.abridged\n\n";

my $abridged_file = $ARGV[0] or die $usage;


=format_input

0       #geneA
1       chr_brkpt_A
2       geneB
3       chr_brkpt_B
4       junction_count
5       spanning_count
6       num_left_contrary_reads
7       num_right_contrary_reads
8       TAF_left
9       TAF_right
10      fusion_annotations
11      TrinGG_Fusion

=cut


main: {

    my @fusion_candidates;


    print STDERR "-parsing $abridged_file\n";
    open (my $fh, $abridged_file) or die "Error, cannot open file $abridged_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my ($geneA, $chr_brkpt_A, 
            $geneB, $chr_brkpt_B, $splice_type,
            $junction_count, $spanning_count, 
            $num_left_contrary_reads, $num_right_contrary_reads,
            $TAF_left, $TAF_right,
            $fusion_annotations, $TrinGG_Fusion) = split(/\t/);

        my $struct = { geneA => $geneA,
                       chr_brkpt_A => $chr_brkpt_A,
                                              
                       geneB => $geneB,
                       chr_brkpt_B => $chr_brkpt_B,

                       splice_type => $splice_type,
                       
                       junction_count => $junction_count,
                       spanning_count => $spanning_count,
                       
                       num_left_contrary_reads => $num_left_contrary_reads,
                       num_right_contrary_reads => $num_right_contrary_reads,

                       TAF_left => $TAF_left,
                       TAF_right => $TAF_right,
                       
                       fusion_annotations => $fusion_annotations,
                       TrinGG_Fusion => $TrinGG_Fusion,
        
                       score => sqrt($junction_count**2 + $spanning_count**2), 
                       
        };
        
        push (@fusion_candidates, $struct);

    }
    close $fh;


    
    ## Filter using homology data
    print STDERR "-parsing blast homology info: $FUSION_ANNOTATOR_LIB/blastn.gene_pairs.gz\n";
    my %blast_pairs = &get_blast_pairs("$FUSION_ANNOTATOR_LIB/blastn.gene_pairs.gz");
    
    print STDERR "-labeling fusions, generating final report.\n";

    @fusion_candidates = &label_likely_artifacts(\@fusion_candidates, \%blast_pairs);


    ## Generate final report
    @fusion_candidates = reverse sort {$a->{score}<=>$b->{score}} @fusion_candidates;

    # header
    print join("\t", "#fusion_name", "score", 
               "geneA", "chr_brkpt_A",
               "geneB", "chr_brkpt_B", "splice_type",
               "junction_count", "spanning_count", 
               "num_left_contrary_reads", "num_right_contrary_reads",
               "TAF_left", "TAF_right",
               "fusion_annots", "TrinityGG_assembled") . "\n";

    foreach my $fusion_candidate (@fusion_candidates) {
        
        my $fusion_name = join("--", $fusion_candidate->{geneA}, $fusion_candidate->{geneB});

        print join("\t", $fusion_name, $fusion_candidate->{score}, 
                   $fusion_candidate->{geneA}, $fusion_candidate->{chr_brkpt_A},
                   $fusion_candidate->{geneB}, $fusion_candidate->{chr_brkpt_B}, 
                   $fusion_candidate->{splice_type},
                   $fusion_candidate->{junction_count}, $fusion_candidate->{spanning_count},
                   $fusion_candidate->{num_left_contrary_reads}, $fusion_candidate->{num_right_contrary_reads}, 
                   $fusion_candidate->{TAF_left}, $fusion_candidate->{TAF_right},
                   $fusion_candidate->{fusion_annotations}, $fusion_candidate->{TrinGG_Fusion}) . "\n";

    }

    print STDERR "-done.\n\n";
    

    exit(0);
}


####
sub label_likely_artifacts {
    my ($fusion_candidates_aref, $blast_pairs_href) = @_;
    
    # sort by score
    my @fusion_candidates = reverse sort {$a->{score}<=>$b->{score}} @$fusion_candidates_aref;
    


    my %seen;
    foreach my $struct (@fusion_candidates) {
        my $geneA = $struct->{geneA};
        my $geneB = $struct->{geneB};

        my $is_likely_artifact = 0;

        if ($seen{$geneA} || $seen{$geneB}) {
            
            if (my $chosenB = $seen{$geneA}) {
                
                if ($blast_pairs_href->{$chosenB}->{$geneB}) {
                    $is_likely_artifact = 1;
                    
                    if ($struct->{fusion_annotations} eq '.') {
                        $struct->{fusion_annotations} = "POTENTIAL_FUSION_ARTIFACT";
                    }
                    else {
                        # append to existing annots
                        $struct->{fusion_annotations} .= ",POTENTIAL_FUSION_ARTIFACT";
                    }
                    $struct->{score} = -1;
                }
            }
                

            if (my $chosenA = $seen{$geneB}) {

                if ($blast_pairs_href->{$chosenA}->{$geneA}) {
                    $is_likely_artifact = 1;
                }
            }
        }

        if ($is_likely_artifact) {
            
            
            
        }
        else {
            $seen{$geneA} = $geneB;
            $seen{$geneB} = $geneA;
        }

    }

    return(@fusion_candidates);
        
}


####
sub get_blast_pairs {
    my ($blast_pairs_file) = @_;

    my %blast_pairs;

    open (my $fh, "gunzip -c $blast_pairs_file | ") or die "Error, cannot open file $blast_pairs_file";
    while (<$fh>) {
        chomp;
        my ($geneA, $geneB, @rest) = split(/\t/);
        
        $blast_pairs{$geneA}->{$geneB} = 1;
        $blast_pairs{$geneB}->{$geneA} = 1;

    }
    close $fh;

    return(%blast_pairs);
}

