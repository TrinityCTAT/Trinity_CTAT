#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Data::Dumper;
use Pipeliner;
use Fastq_reader;

my $usage = <<__EOUSAGE__;

###############################################################################################
#
#  --left_fq <string>             left fastq filename
#  
#  --right_fq <string>            right fastq filename
#
#  --cdna_fa <string>             reference cdna fasta filename
#
#  --fusion_summary <string>      fusion summary file (pre-filtering)
#
###############################################################################################


__EOUSAGE__

    ;

my $left_fq;
my $right_fq;
my $cdna_fa;
my $fusion_summary;
my $help_flag;

&GetOptions ( 'h' => \$help_flag,
              
              'left_fq=s' => \$left_fq,
              'right_fq=s' => \$right_fq,
              'cdna_fa=s' => \$cdna_fa,
              'fusion_summary=s' => \$fusion_summary,
    );


if ($help_flag) {
    die $usage;
}

unless ($left_fq && $right_fq && $cdna_fa && $fusion_summary) {
    die $usage;
}


main: {


    print STDERR "-parsing junction and spanning reads info.\n";
    my %junction_reads;
    my %spanning_frags;
    &parse_junction_and_spanning_reads($fusion_summary, \%junction_reads, \%spanning_frags);


    {
        ##
        ## Examine/filter the junction reads
        ##

        print STDERR "-extracting junction reads from fq file\n";
        
        my $junc_reads_fq = "tmp.junc_reads.fq";
        my $junc_reads_chkpt = "$junc_reads_fq.ok";
        ## examine junction reads
        unless (-e $junc_reads_chkpt) {
            
            ## get the reads
            &extract_junction_reads(\%junction_reads, $left_fq, $right_fq, $junc_reads_fq);
            
            
            ## align the reads to the cdna fasta file
            
            
            
            &process_cmd("touch $junc_reads_chkpt");
            
        }
    }

    {
        ##
        ## Examine/filter the spanning frags
        ##

        print STDERR "-extracting spanning frags from fq files\n";
        
        my $span_reads_left_fq = "tmp.span_reads.left.fq";
        my $span_reads_right_fq = "tmp.span_reads.right.fq";
        
        my $span_reads_chkpt = "tmp.span_reads.ok";
        unless (-e $span_reads_chkpt) {

            &extract_spanning_reads(\%spanning_frags, $left_fq, $span_reads_left_fq);
            &extract_spanning_reads(\%spanning_frags, $right_fq, $span_reads_right_fq);
                        

            &process_cmd("touch $span_reads_chkpt");
        }
     
        ## align reads, identify those that actually align as proper pairs.
        
   
    }

    exit(0);

}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        confess "Error, cmd: $cmd died with ret $ret";
    }
    return;
}


####
sub extract_junction_reads {
    my ($junction_reads_href, $left_fq, $right_fq, $junc_reads_fq) = @_;

    my %junc_reads = %$junction_reads_href;

    open (my $ofh, ">$junc_reads_fq") or confess "Error, cannot write to file $junc_reads_fq";

    foreach my $fq_file ($left_fq, $right_fq) {
        
        my $fastq_reader = new Fastq_reader($fq_file);
        while (my $fq_obj = $fastq_reader->next()) {
            my $full_read_name = $fq_obj->get_full_read_name();
            if (exists $junc_reads{$full_read_name}) {
                my $fq_record = $fq_obj->get_fastq_record();
                print $ofh $fq_record;
                delete $junc_reads{$full_read_name};
            }
        }
    }
    close $ofh;
    
    if (%junc_reads) {
        confess "Error, didn't retrieve all junction reads.  Missing: " . Dumper(\%junc_reads);
    }

    return;
}

####
sub extract_spanning_reads {
    my ($spanning_frags_href, $fastq_input, $fastq_output) = @_;

    my %frags = %$spanning_frags_href;
    
    my $fastq_reader = new Fastq_reader($fastq_input);
    open (my $ofh, ">$fastq_output") or die "error, cannot write to file $fastq_output";
    while (my $fq_obj = $fastq_reader->next()) {
        my $core_read_name = $fq_obj->get_core_read_name();
        if (exists $frags{$core_read_name}) {
            my $fq_record = $fq_obj->get_fastq_record();
            print $ofh $fq_record;
            delete($frags{$core_read_name});
        }
    }
    close $ofh;

    return;
}


####
sub parse_junction_and_spanning_reads {
    my ($fusion_summary, $junction_reads_href, $spanning_frags_href) = @_;

    open (my $fh, $fusion_summary) or die "Error, cannot open file $fusion_summary";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);
        my $junc_reads_list = $x[9];
        my $spanning_frag_list = $x[10];

        foreach my $junc_read (split(/,/, $junc_reads_list)) {
            $junction_reads_href->{$junc_read}++;
        }
        foreach my $spanning_frag (split(/,/, $spanning_frag_list)) {
            $spanning_frags_href->{$spanning_frag}++;
        }
    }
    close $fh;

    return;
}
