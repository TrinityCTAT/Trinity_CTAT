#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use FusionAnnotator;

my $usage = "\n\n\tusage: $0 Chimeric.out.junction.genes Chimeric.out.sam.spanning.genes\n\n";

my $junction_genes = $ARGV[0] or die $usage;
my $spanning_genes = $ARGV[1] or die $usage;

main: {

    my %fusion_to_counts;

    my %fusion_to_coords_counts;

    my %fusion_read_names;

    {
        open (my $fh, $junction_genes) or die "Error, cannot open file $junction_genes";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $frag_name = $x[9];
            my $fusion = $x[14];
            my $fusion_coords_text = join("\t", @x[15..$#x]);
            
            if ($fusion =~ /unknown/) { next; }
            
            $fusion_to_counts{$fusion}->{junction}++;
            
            $fusion_to_coords_counts{$fusion}->{$fusion_coords_text}++;
        
            $fusion_read_names{$frag_name}++;
            
        }
        close $fh;
    }

    {
        open (my $fh, $spanning_genes) or die "Error, cannot open file $spanning_genes";
        while (<$fh>) {
            chomp;
            my ($frag_name, $fusion) = split(/\t/);

            if ($fusion_read_names{$frag_name}) { next; } # no over-counting

            my ($genesA, $genesB) = split(/--/, $fusion);
            
            my %seen;
            foreach my $geneA (split(/,/, $genesA)) {
                foreach my $geneB (split(/,/, $genesB)) {
                    
                    my $fusion_A = "$geneA--$geneB";
                    my $fusion_B = "$geneB--$geneA";
                    
                    if ($seen{$fusion_A}) { next; } # no over-counting
                    $seen{$fusion_A} = 1;

                    if (exists $fusion_to_counts{$fusion_A}) {
                        $fusion_to_counts{$fusion_A}->{spanning}++;
                    }
                    if (exists $fusion_to_counts{$fusion_B}) {
                        $fusion_to_counts{$fusion_B}->{spanning}++;
                    }
                }
            }
        }
        
        close $fh;
        
    }
    
    

    foreach my $fusion (keys %fusion_to_counts) {


        my ($geneA, $geneB) = split(/--/, $fusion);
        
        ## no self-fusions
        if ($geneA eq $geneB) { next; }


        my @annotations = &FusionAnnotator::get_annotations($geneA, $geneB);

        my $junction_count = $fusion_to_counts{$fusion}->{junction};
        my $spanning_count = $fusion_to_counts{$fusion}->{spanning} || 0;
        
        my $fusion_coords_info_href = $fusion_to_coords_counts{$fusion};
        
        # get the coordinate junction set with the highest read support
        my @coords_info = reverse sort {$fusion_coords_info_href->{$a} <=> $fusion_coords_info_href->{$b}} keys %$fusion_coords_info_href;

        print join("\t", $fusion, $junction_count, $spanning_count, $coords_info[0],
                   join(",", @annotations)
            ) . "\n";
    }
    
    exit(0);
}
            
