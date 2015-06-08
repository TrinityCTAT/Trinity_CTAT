#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use __GLOBALS__;
use FusionAnnotator;
use Data::Dumper;

my $MAX_PCT_J_EV_READS_FILTERED = 90;

my $MIN_SCORE_PCT_FUSION_TO_PARALOG_FUSION = 90; # if there exists fusion A-B but also a predicted A-C where B,C are paralogs and A-C scores <= A-B, store A-C as an alternate under A-B.

my $min_novel_junction_support = 10; # same as we use in STAR-Fusion


my $usage = "\n\tusage: $0 fi_test.fusion_preds.coalesced.summary.wTrinityGG.abridged\n\n";

my $abridged_file = $ARGV[0] or die $usage;


=format_input

0       #geneA
1       chr_brkpt_A
2       geneB
3       chr_brkpt_B
4       splice_type
5       junction_count
6       spanning_count
7       num_left_contrary_reads
8       num_right_contrary_reads
9       TAF_left
10       TAF_right
11      fusion_annotations
12      TrinGG_Fusion

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

        unless ($TrinGG_Fusion) {
            $TrinGG_Fusion = "."; # adding placeholder
        }

        my $struct = { 

            fusion_name => "$geneA--$geneB",
            
            geneA => $geneA,
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
            
            alternates => [],
            
        };
        
        push (@fusion_candidates, $struct);

    }
    close $fh;


    @fusion_candidates = &label_likely_artifacts(\@fusion_candidates);
    
    
    ## Generate final report
    @fusion_candidates = reverse sort {$a->{score}<=>$b->{score}} @fusion_candidates;

    # header
    print join("\t", "#fusion_name", "score", 
               "geneA", "chr_brkpt_A",
               "geneB", "chr_brkpt_B", "splice_type",
               "junction_count", "spanning_count", 
               "num_left_contrary_reads", "num_right_contrary_reads",
               "TAF_left", "TAF_right",
               "fusion_annots", "TrinityGG_assembled", "Alternates") . "\n";

    foreach my $fusion_candidate (@fusion_candidates) {
        
        my $fusion_name = join("--", $fusion_candidate->{geneA}, $fusion_candidate->{geneB});

        my $alternate_fusions = join(",", @{$fusion_candidate->{alternates}}) || ".";

        print join("\t", $fusion_name, sprintf("%.2f", $fusion_candidate->{score}), 
                   $fusion_candidate->{geneA}, $fusion_candidate->{chr_brkpt_A},
                   $fusion_candidate->{geneB}, $fusion_candidate->{chr_brkpt_B}, 
                   $fusion_candidate->{splice_type},
                   $fusion_candidate->{junction_count}, $fusion_candidate->{spanning_count},
                   $fusion_candidate->{num_left_contrary_reads}, $fusion_candidate->{num_right_contrary_reads}, 
                   $fusion_candidate->{TAF_left}, $fusion_candidate->{TAF_right},
                   $fusion_candidate->{fusion_annotations}, $fusion_candidate->{TrinGG_Fusion},
                   $alternate_fusions,
                   
            ) . "\n";

    }
    
    print STDERR "-done.\n\n";
    

    exit(0);
}


####
sub label_likely_artifacts {
    my ($fusion_candidates_aref) = @_;
    
    # sort by score
    my @fusion_candidates = reverse sort {$a->{score}<=>$b->{score}} @$fusion_candidates_aref;
    

    #print Dumper(\@fusion_candidates);
    

    my %AtoB;
    my %BtoA;
    
    my %selected;
    
    my %fusion_name_to_fusion_struct;

    foreach my $struct (@fusion_candidates) {
        my $geneA = $struct->{geneA};
        my $geneB = $struct->{geneB};

        my $fusion_name = "$geneA--$geneB";
        unless (exists $fusion_name_to_fusion_struct{$fusion_name}) {
            $fusion_name_to_fusion_struct{$fusion_name} = $struct;
        }
        
        my $score = $struct->{score};
        my $fusion_annot = $struct->{fusion_annotations};

        my $is_likely_artifact = 0;
        
        my $fusion_splice_type = $struct->{splice_type};
        my $num_junction_reads = $struct->{junction_count};

        if ($fusion_splice_type eq "INCL_NON_REF_SPLICE" && $num_junction_reads < $min_novel_junction_support) {
            
            $is_likely_artifact = 1;
            
        }
        elsif($fusion_annot =~ /PctFiltJ\[([^\]]+)/g) {
            # PctFiltJ[$pct_filtered_junction],PctFiltS[$pct_filtered_spanning]
            my $pct_filt = $1;
            if ($pct_filt > $MAX_PCT_J_EV_READS_FILTERED) {
                $is_likely_artifact = 1;
                last;
            }
        }
        
        if ( (! $is_likely_artifact) && ($AtoB{$geneA} || $BtoA{$geneB}) ) {
            
            if (my $chosenB_href = $AtoB{$geneA}) {
                
                foreach my $chosenB (keys %$chosenB_href) {
                    
                    if ( grep { /BLAST|CLUST/ } &FusionAnnotator::get_annotations($geneB, $chosenB) ) {

                        $is_likely_artifact = 1;
                    
                        my $chosen_struct = $fusion_name_to_fusion_struct{"$geneA--$chosenB"} or confess "Error, no struct stored for fusion: $geneA--$chosenB";
                        
                        if ($score >= $MIN_SCORE_PCT_FUSION_TO_PARALOG_FUSION/100 * $chosen_struct->{score}) {
                            
                            push (@{$chosen_struct->{alternates}}, $fusion_name);
                        
                        }
                        
                    }
                    
                }
            }
            
            
            if (my $chosenA_href = $BtoA{$geneB}) {
                
                foreach my $chosenA (keys %$chosenA_href) {
                   
                    if (grep { /BLAST|CLUST/ } &FusionAnnotator::get_annotations($geneA, $chosenA) ) {
                        
                        $is_likely_artifact = 1;
                        
                        my $chosen_struct = $fusion_name_to_fusion_struct{"$chosenA--$geneB"} or confess "Error, no struct stored for fusion: $chosenA--$geneB";
                        
                        if ($score >= $MIN_SCORE_PCT_FUSION_TO_PARALOG_FUSION/100 * $chosen_struct->{score}) {
                            
                            push (@{$chosen_struct->{alternates}}, $fusion_name);
                            
                        }                    
                        
                    }
                }
            }
        }
        
        
        if ($is_likely_artifact) {
            if ($struct->{fusion_annotations} eq '.') {
                $struct->{fusion_annotations} = "POTENTIAL_FUSION_ARTIFACT";
            }
            else {
                # append to existing annots
                $struct->{fusion_annotations} .= ",POTENTIAL_FUSION_ARTIFACT";
            }
            $struct->{score} = -1;
        }
                
        else {
            $AtoB{$geneA}->{$geneB} = 1;
            $BtoA{$geneB}->{$geneA} = 1;
            
        }
        
    }

    return(@fusion_candidates);
        
}
