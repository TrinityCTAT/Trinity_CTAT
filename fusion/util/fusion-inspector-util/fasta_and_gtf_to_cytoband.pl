#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 genome.fasta annotations.gtf\n\n";

my $genome_fasta = $ARGV[0] or die $usage;
my $annotations = $ARGV[1] or die $usage;

main: {


    my %contig_to_genes = &parse_gtf($annotations);
    
    print Dumper(\%contig_to_genes);

    

    my $fasta_reader = new Fasta_reader($genome_fasta);
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        my $seq_len = length($sequence);

        if (! exists $contig_to_genes{$acc}) {
            die "Error, no genes on $acc";
        }

        my @genes = @{$contig_to_genes{$acc}};
                

        @genes = sort {$a->{lend}<=>$b->{lend}} @genes;

        my $prev_lend = 1;
        foreach my $gene (@genes) {
            my $gene_id = $gene->{gene_id};
            my $lend = $gene->{lend};
            my $rend = $gene->{rend};
            
            print join("\t", $acc, $prev_lend, $lend, "spacer_$prev_lend", "gneg") . "\n";
            print join("\t", $acc, $lend, $rend, $gene_id, "gpos100") . "\n";

            $prev_lend = $rend;
        }
        
        # get last one.
        print join("\t", $acc, $prev_lend, $seq_len, "spacer_$prev_lend", "gneg") . "\n";
    }


    exit(0);



}


####
sub parse_gtf {
    my ($annotations) = @_;

    my %contig_to_gene_coords;
    
    open (my $fh, $annotations) or die "Error, cannot open file $annotations";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $contig = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        
        my $info = $x[8];

        $info =~ /gene_id \"([^\"]+)/ or die "Error, cannot parse gene_id info from $info";
        my $gene_id = $1;

        if ($info =~ /gene_name \"([^\"]+)/) {
            # use gene name instead
            $gene_id = $1;
        }

        $gene_id = "$contig|$gene_id";
        
        push (@{$contig_to_gene_coords{$gene_id}}, $lend, $rend);

    }
    close $fh;

    my %contig_to_genes;

    foreach my $gene_id (keys %contig_to_gene_coords) {
        
        my ($contig, $core_gene_id) = split(/\|/, $gene_id);

        my @coords = sort {$a<=>$b} @{$contig_to_gene_coords{$gene_id}};
        my $gene_lend = shift @coords;
        my $gene_rend = pop @coords;

        
        my $struct = { gene_id => $gene_id,
                       lend => $gene_lend,
                       rend => $gene_rend,
        };
        
        push (@{$contig_to_genes{$contig}}, $struct);
    }

    return(%contig_to_genes);
}


