#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw (min max);
use POSIX;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_retriever;
use Fasta_reader;

my $usage = "usage: $0  gmap.map.gff3.chims_described gmap.map.gff3.chims_described.fasta reference_cdna.fasta KMER_SIZE ANCHOR_SIZE\n\n";

my $fusion_preds_file = $ARGV[0] or die $usage;
my $fusion_trans_file = $ARGV[1] or die $usage;
my $ref_cdna_fasta = $ARGV[2] or die $usage;
my $KMER_SIZE = $ARGV[3] or die $usage;
my $ANCHOR = $ARGV[4] or die $usage;


main: {
    
    my $fasta_reader = new Fasta_reader($fusion_trans_file);
    my %fusion_trans = $fasta_reader->retrieve_all_seqs_hash();
    
    my %gene_to_trans_list = &parse_trans_to_gene_list($ref_cdna_fasta);
    
    my $fasta_retriever = new Fasta_retriever($ref_cdna_fasta);
    
    open (my $fh, $fusion_preds_file) or die "Error, cannot open file $fusion_preds_file";
    open (my $ofh, ">$fusion_preds_file.kmer_confound_info") or die "Error, cannot write to file $fusion_preds_file.kmer_confound_info";

    while (<$fh>) {
        if (/^\#/) { next; }

        chomp;
        my $line = $_;
        my @x = split(/\t/);
        
        my $trans_acc = $x[0];
        my $fusion_info = $x[3];
        
        my @fusion_pts = split(/;/, $fusion_info);
        
        my $geneA = $fusion_pts[0];
        my $trans_brk_left = $fusion_pts[2];
        
        my $geneB = $fusion_pts[4];
        my $trans_brk_right = $fusion_pts[6];


        my $left_gene_seqs = &get_gene_seqs($geneA, \%gene_to_trans_list, $fasta_retriever);
        my $right_gene_seqs = &get_gene_seqs($geneB, \%gene_to_trans_list, $fasta_retriever);
        
                
        my %gene_shared_kmers = &get_sequence_shared_kmers($left_gene_seqs, $right_gene_seqs);
        
        ## see if the shared kmers occur near the breakpoint of the candidate fusion transcript
        my $fusion_trans_seq = $fusion_trans{$trans_acc};
        my %fusion_breakpoint_kmers = &get_kmers($fusion_trans_seq, $trans_brk_left-$ANCHOR, $trans_brk_right+$ANCHOR);
        
        my %shared_breakpoint_kmers = &subset_kmers(\%gene_shared_kmers, \%fusion_breakpoint_kmers);
        
        if (%shared_breakpoint_kmers) {
            
            ## filter and annotate
            my $kmer_text = "";
            foreach my $kmer (keys %shared_breakpoint_kmers) {
                my @kmer_pos_list = @{$shared_breakpoint_kmers{$kmer}};
                $kmer_text .= "\t$kmer(" . join(",", @kmer_pos_list) . ")";
            }
            print $ofh "#$line" . $kmer_text . "\n";
        }
        else {
            print $ofh $line . "\n";
            print $line . "\n"; # unfiltered goes to stdout.
        }
        
    }
    close $fh;
    close $ofh;


    exit(0);
    
}

####
sub get_sequence_shared_kmers {
    my ($seqA, $seqB) = @_;
    
    my %kmersA = &get_kmers($seqA);

    my %kmersB = &get_kmers($seqB);
    
    my %joint_kmers;
    foreach my $kmer (keys %kmersA) {
        if (exists $kmersB{$kmer}) {
            my $A_pos_list = $kmersA{$kmer};
            my $B_pos_list = $kmersB{$kmer};
            
            $joint_kmers{$kmer} = [$A_pos_list, $B_pos_list];
        }
    }

    return(%joint_kmers);
}


####
sub get_kmers {
    my ($seq, $beg_pos, $end_pos)= @_;

    my %kmers_to_pos_list;
        
    unless (defined $beg_pos) {
        $beg_pos = 0;
    }
    unless (defined $end_pos) {
        $end_pos = length($seq);
    }
    if ($end_pos > length($seq)) {
        $end_pos = length($seq);
    }
    
    unless ($end_pos >= $beg_pos) {
        die "Error, beginning and end positions are out of order: beg($beg_pos) and end($end_pos)";
    }
    
    # print "WINDOW ($centerpt): $beg_pos to $end_pos\n";
    
    for (my $i = $beg_pos; $i <= $end_pos - $KMER_SIZE; $i++) {
        
        my $kmer = uc substr($seq, $i, $KMER_SIZE);
        
        push (@{$kmers_to_pos_list{$kmer}}, $i);
    }

    return(%kmers_to_pos_list);
}


####
sub parse_trans_to_gene_list {
    my ($ref_cdna_fasta) = @_;

    my %gene_to_trans;

    open (my $fh, $ref_cdna_fasta) or die "Error, cannot open file $ref_cdna_fasta";
    while (<$fh>) {
        chomp;
        if (/^>/) {
            s/>//;
            my ($trans_id, $gene_id, $gene_name) = split(/\s+/);
            push (@{$gene_to_trans{$gene_id}}, $trans_id);
            if ($gene_name && $gene_name ne $gene_id) {
                push (@{$gene_to_trans{$gene_name}}, $trans_id);
            }
        }
    }

    close $fh;


    return(%gene_to_trans);
}


####
sub get_gene_seqs {
    my ($gene, $gene_to_trans_list_href, $fasta_retriever) = @_;
    
    unless (exists $gene_to_trans_list_href->{$gene}) {
        die "Error, no list of transcript ids for gene: $gene";
    }
    my @trans_ids = @{$gene_to_trans_list_href->{$gene}};

    my @trans_seqs;
    foreach my $trans_id (@trans_ids) {
        my $seq = $fasta_retriever->get_seq($trans_id) or die "Error, no sequence found for trancript $trans_id";
        push (@trans_seqs, $seq);
    }

    my $concat_gene_seq = join("X", @trans_seqs);
    
    return($concat_gene_seq);
}
    
####
sub subset_kmers {
    my ($kmer_hashref, $kmer_subset_hashref) = @_;

    my %kmers_want;

    foreach my $kmer (keys %$kmer_subset_hashref) {
        if (exists $kmer_hashref->{$kmer}) {
            $kmers_want{$kmer} = $kmer_subset_hashref->{$kmer};
        }
    }

    return(%kmers_want);
}

