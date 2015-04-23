#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use __GLOBALS__;
use FusionAnnotator;

my $usage = "\n\tusage: $0 junction_info_A.txt,[junction_info_B.txt,...] spanning_info_A.txt,[spanning_info_B.txt,...]\n\n";

my $junction_info_file_list = $ARGV[0] or die $usage;
my $spanning_info_file_list = $ARGV[1] or die $usage;

my $PSEUDOCOUNT = 1;

main: {

    my %fusion_info;

    foreach my $junction_info_file (split(/,/, $junction_info_file_list)) {

        &parse_info_file(\%fusion_info, "junction", $junction_info_file);
    }
    
    foreach my $spanning_info_file (split(/,/, $spanning_info_file_list)) {

        &parse_info_file(\%fusion_info, "spanning", $spanning_info_file);

    }

    # generate coalesced view:
    print join("\t", "#geneA", "local_brkpt_A", "chr_brkpt_A", 
               "geneB", "local_brkpt_B", "chr_brkpt_B", "splice_type",
               "junction_count", "spanning_count", "junction_reads", "spanning_reads", 
               "num_left_contrary_reads", "left_contrary_reads",
               "num_right_contrary_reads", "right_contrary_reads",
               "TAF_left", "TAF_right", "fusion_annotations") . "\n";

    foreach my $fusion (keys %fusion_info) {
        my @junction_reads = keys %{$fusion_info{$fusion}->{'junction'}};

        my @spanning_reads;
        if (exists $fusion_info{$fusion}->{'spanning'}) {
            @spanning_reads = keys %{$fusion_info{$fusion}->{'spanning'}};
        
            @spanning_reads = &remove_residual_junction_from_spanning_reads(\@junction_reads, \@spanning_reads);

        }

        my $num_junction_reads = scalar(@junction_reads);
        my $num_spanning_reads = scalar(@spanning_reads);

        my @left_contrary_reads;
        if (exists $fusion_info{$fusion}->{'spanning-left_contrary_support'}) {
            @left_contrary_reads = keys %{$fusion_info{$fusion}->{'spanning-left_contrary_support'}};
        }
        my $num_left_contrary_reads = scalar(@left_contrary_reads);
        if ($num_left_contrary_reads == 0) {
            push (@left_contrary_reads, ".");
        }
        

        my @right_contrary_reads;
        if (exists $fusion_info{$fusion}->{'spanning-right_contrary_support'}) {
            @right_contrary_reads = keys %{$fusion_info{$fusion}->{'spanning-right_contrary_support'}};
        }
        my $num_right_contrary_reads = scalar(@right_contrary_reads);
        if ($num_right_contrary_reads == 0) {
            push (@right_contrary_reads, ".");
        }
        
        my $TAF_left = ($num_junction_reads + $num_spanning_reads + $PSEUDOCOUNT) / ($num_left_contrary_reads + $PSEUDOCOUNT);
        $TAF_left = sprintf("%.2f", $TAF_left);
        
        my $TAF_right = ($num_junction_reads + $num_spanning_reads + $PSEUDOCOUNT) / ($num_right_contrary_reads + $PSEUDOCOUNT);
        $TAF_right = sprintf("%.2f", $TAF_right);
        
        my ($geneA, $local_brkpt_A, $chr_brkpt_A, $geneB, $local_brkpt_B, $chr_brkpt_B) = split(/\t/, $fusion);
        
        
        my @annots = &FusionAnnotator::get_annotations($geneA, $geneB);
        unless (@annots) {
            @annots = ("."); # put in a placeholder
        }


        print join("\t", $fusion, $num_junction_reads, $num_spanning_reads,
                   join(",", @junction_reads),
                   join(",", @spanning_reads),
                   $num_left_contrary_reads,
                   join(",", @left_contrary_reads),
                   $num_right_contrary_reads,
                   join(",", @right_contrary_reads),
                   $TAF_left,
                   $TAF_right,
                   join(",", @annots),
            ) 
            
            . "\n";
    }
    


    exit(0);
}


####
sub parse_info_file {
    my ($fusion_info_href, $fusion_read_type, $file) = @_;

    open (my $fh, $file) or die "Error, cannot open file $file";
    while (<$fh>) {
        chomp;
        my ($geneA, $coordA, $orig_coordA, $geneB, $coordB, $orig_coordB, $splice_info, $count, $read_list, @rest) = split(/\t/);
        
        my $fusion_token = join("\t", $geneA, $coordA, $orig_coordA, $geneB, $coordB, $orig_coordB, $splice_info);
        
        foreach my $read (split(/,/, $read_list)) {

            $fusion_info_href->{$fusion_token}->{$fusion_read_type}->{$read}++;
        }


        if (@rest) {  # in spanning file
            
            my ($left_contrary_support_count, $left_contrary_support_reads,
                $right_contrary_support_count, $right_contrary_support_reads) = @rest;
            
            if ($left_contrary_support_count > 0) {
                foreach my $read (split(/,/, $left_contrary_support_reads)) {
                    
                    $fusion_info_href->{$fusion_token}->{"$fusion_read_type-left_contrary_support"}->{$read} = 1;
                }
            }
            
            if ($right_contrary_support_count > 0) {
                foreach my $read (split(/,/, $right_contrary_support_reads)) {
                    
                    $fusion_info_href->{$fusion_token}->{"$fusion_read_type-right_contrary_support"}->{$read} = 1;
                    
                }
            }
        }
    }
    close $fh;

    
    return;
}


####
sub remove_residual_junction_from_spanning_reads {
    my ($junction_reads_aref, $spanning_reads_aref) = @_;

    # junction takes precendence over spanning.
    
    my %junction;
    foreach my $junction_read (@$junction_reads_aref) {
        
        my $read_to_discard = $junction_read;
        $read_to_discard =~ s|/[12]$||;
        
        $junction{$read_to_discard}++;
    }

    my @spanning_retain;
    
    foreach my $spanning_read (@$spanning_reads_aref) {
        if (exists $junction{$spanning_read}) {
            print STDERR "-discarding spanning read: $spanning_read since found among merged junction reads.\n";
        }
        else {
            push (@spanning_retain, $spanning_read);
        }
    }

    return(@spanning_retain);
}


