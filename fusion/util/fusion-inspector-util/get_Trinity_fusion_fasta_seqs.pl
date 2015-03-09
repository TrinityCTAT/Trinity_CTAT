#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;


my $usage = "\n\tusage: $0 Trinty.fasta Trinity.fusions.gff3\n\n";

my $trin_fasta = $ARGV[0] or die $usage;
my $trin_fusions_gff3 = $ARGV[1] or die $usage;


## read the Trinity accessions and breakpoint info from the header comments of the gff3


main: {

    my %trin_acc_to_breakpoint_info;
    {
        open (my $fh, $trin_fusions_gff3) or die "Error, cannot open file $trin_fasta";
        while (<$fh>) {
            chomp;
            if (/^\#TrinityFusionTranscript/) {
                chomp;
                my @x = split(/\t/);
                my ($token, $trin_acc, $contig_brkpt) = @x;

                $trin_acc_to_breakpoint_info{$trin_acc} = $contig_brkpt;
            }
        }
        close $fh;
    }

    my $fasta_reader = new Fasta_reader($trin_fasta);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        if (my $contig_brkpt = $trin_acc_to_breakpoint_info{$acc}) {
            
            my $sequence = $seq_obj->get_sequence();

            print ">$acc $contig_brkpt\n$sequence\n";
        }
    }
    

    exit(0);
}


            
