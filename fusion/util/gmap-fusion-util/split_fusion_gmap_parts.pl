#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;


my $usage = "\n\n\tusage: $0 gmap.map.gff3.chims_described gmap.map.gff3.chims_described.fasta\n\n";

my $chims_described_file = $ARGV[0] or die $usage;
my $chims_fasta_file = $ARGV[1] or die $usage;


## configuration:
my $FUSION_ANNOTATOR_LIB = $ENV{FUSION_ANNOTATOR_LIB} or die "Error, must specify env var FUSION_ANNOTATOR_LIB";
my $GENOME = "${FUSION_ANNOTATOR_LIB}/Hg19.fa";
my $GMAP_DB_DIR = "${FUSION_ANNOTATOR_LIB}";
my $GMAP_DB_NAME = "Hg19.fa.gmap";


main: {

    my %transcript_to_breakpoint = &parse_chimera_preds($chims_described_file);

    my $fasta_reader = new Fasta_reader($chims_fasta_file);
    my %trans_seqs = $fasta_reader->retrieve_all_seqs_hash();

    my $chim_frag_file = "$chims_fasta_file.split.fa";
    open (my $ofh, ">$chim_frag_file");

    foreach my $trans (keys %transcript_to_breakpoint) {
        
        my $sequence = $trans_seqs{$trans};
        
        my $breakpoint_range = $transcript_to_breakpoint{$trans};
        my ($brk_left, $brk_right) = sort {$a<=>$b} split(/-/, $breakpoint_range);

        my $seq_range_left = substr($sequence, 0, $brk_left);
        my $seq_range_right = substr($sequence, $brk_right);

        print $ofh ">$trans" . "____left\n"
            . "$seq_range_left\n"
            . ">$trans" . "____right\n"
            . "$seq_range_right\n";
        
    }
    close $ofh;

    ## run GMAP, capture all top hits within reason.
    my $gmap_output_file = "$chim_frag_file.gmap.gff3";
    my $cmd = "gmap -D $GMAP_DB_DIR -d $GMAP_DB_NAME $chim_frag_file -f 3 -n 10 -t 4  > $gmap_output_file";
    
    my $pipeliner = new Pipeliner(-verbose => 2);
    $pipeliner->add_commands(new Command($cmd, "gmap_output_file.ok"));

    $pipeliner->run();
    
    
    exit(0);
}


####
sub parse_chimera_preds {
    my ($chims_described_file) = @_;

    my %trans_to_brk;

    open (my $fh, $chims_described_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $trans_acc = $x[0];
        my $info = $x[2];
        my @pts = split(/;/, $info);
        my $brk_left = $pts[2];
        my $brk_right = $pts[6];

        $trans_to_brk{$trans_acc} = "$brk_left-$brk_right";

    }

    close $fh;

    return(%trans_to_brk);
}
        
