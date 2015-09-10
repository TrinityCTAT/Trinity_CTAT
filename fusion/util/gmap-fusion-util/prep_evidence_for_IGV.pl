#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Process_cmd;

my $usage = "\n\n\tusage: $0 gmap.map.gff3.chims_described.fasta gmap.map.gff3.chims_described.w_read_support.J1.JS2.NJ3 gmap.map.gff3.chims_described.fasta.bowtie2.bam\n\n";


my $chims_fasta = $ARGV[0] or die $usage;
my $fusion_preds = $ARGV[1] or die $usage;
my $bam_file = $ARGV[2] or die $usage;

main: {

    unless (-e "$chims_fasta.fai") {
        &process_cmd("samtools faidx $chims_fasta");
    }

    my %junction_reads;
    my %spanning_frags;
    my %trans_breakpoint;
    &parse_fusion_evidence($fusion_preds, \%junction_reads, \%spanning_frags, \%trans_breakpoint);

    &write_fusion_evidence_bams($chims_fasta, $bam_file, \%junction_reads, \%spanning_frags);
    
    

}


####
sub write_fusion_evidence_bams {
    my ($chims_fasta, $bam_file, $junction_reads_href, $spanning_frags_href) = @_;

    my $junction_frags_sam_file = "junction_reads.sam";
    my $spanning_frags_sam_file = "spanning_frags.sam";

    open (my $ofh_junc, ">$junction_frags_sam_file") or die $!;
    open (my $ofh_span, ">$spanning_frags_sam_file") or die $!;

    my $sam_reader = new SAM_reader($bam_file);
    while (my $sam_entry = $sam_reader->get_next()) {
        my $trans_acc = $sam_entry->get_scaffold_name();
        my $frag_name = $sam_entry->get_read_name();
        
        my $token = join("$;", $trans_acc, $frag_name);

        if (exists $junction_reads_href->{$token}) {
            print $ofh_junc $sam_entry->get_original_line() . "\n";
        }
        elsif (exists $spanning_frags_href->{$token}) {
            print $ofh_span $sam_entry->get_original_line() . "\n";
        }
    }

    close $ofh_junc;
    close $ofh_span;

    
    ## convert to coord-sorted bam files.




    return;

}




####
sub parse_fusion_evidence {
    my ($fusion_preds, $junction_reads_href, $spanning_frags_href, $trans_breakpoint_href) = @_;

    open (my $fh, $fusion_preds) or die "Error, cannot open file $fusion_preds";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $trans_acc = $x[0];
        
        my $fusion_info = $x[3];
        my @pts = split(/;/, $fusion_info);
        my $brkpt = join("-", $pts[2], $pts[7]);
        $trans_breakpoint_href->{$trans_acc} = $brkpt;

        my $junction_read_list = $x[6];
        my $spanning_frag_list = $x[7];
        
        

        foreach my $junction_read (split(/,/, $junction_read_list)) {
            my $token = join("$;", $trans_acc, $junction_read);

            $junction_reads_href->{$token} = 1;
        }

        foreach my $spanning_frag (split(/,/, $spanning_frag_list) ) {
            my $token = join("$;", $trans_acc, $spanning_frag);
            
            $spanning_frags_href->{$token} = 1;
        }

    }
    close $fh;

    
    return;
}
