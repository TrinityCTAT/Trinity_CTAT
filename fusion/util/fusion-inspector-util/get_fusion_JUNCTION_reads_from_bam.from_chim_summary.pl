#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Carp;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 genePairContig.gtf read_alignments.bam\n\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $bam_file = $ARGV[1] or die $usage;

main: {
    
    my %exon_bounds;
    my %original_exon_coord_mapping;
    my %scaffold_to_gene_structs = &parse_gtf_file($gtf_file, \%exon_bounds, \%original_exon_coord_mapping);
    
    print STDERR Dumper(\%scaffold_to_gene_structs);
    
    my $prev_scaff = "";

    my ($left_gene, $right_gene, $gene_bound_left, $gene_bound_right);
    

    my %fusion_junctions;

    ## find the reads that matter:
    my $sam_reader = new SAM_reader($bam_file);
    while (my $sam_entry = $sam_reader->get_next()) {
        
        my $read_name = $sam_entry->reconstruct_full_read_name();

        my $scaffold = $sam_entry->get_scaffold_name();
        if ($scaffold ne $prev_scaff) {
            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            if (scalar @gene_structs != 2) {
                die "Error, didn't find only 2 genes in the gtf file: " . Dumper(\@gene_structs);
            }
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;
            
            $left_gene = $gene_structs[0];
            $right_gene = $gene_structs[1];
            
            $gene_bound_left = $left_gene->{rend};
            $gene_bound_right = $right_gene->{lend};
            
            
            if ($gene_bound_left > $gene_bound_right) {
                die "Error, bounds out of order: $gene_bound_left ! <  $gene_bound_right";
            }
            
            
            $prev_scaff = $scaffold;
        }
        
        my ($span_lend, $span_rend) = $sam_entry->get_genome_span();
        if ($span_lend < $gene_bound_left && $span_rend > $gene_bound_right) {
            
            my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
            
            my @genes_matched = &hits_exon_bound($genome_coords_aref, $exon_bounds{$scaffold});

            if (scalar @genes_matched == 2) {
                
                print $sam_entry->get_original_line() . "\n";
            
                my $junction_coord_token = &get_junction_coord_token($scaffold, \%original_exon_coord_mapping, \@genes_matched);

                push (@{$fusion_junctions{$junction_coord_token}}, $read_name);
            }
        }

    }

    {
        
        print STDERR "-writing fusion junction support info.\n";
        
        ## output the junction info summary:
        my $junction_coords_info_file = "$bam_file.fusion_junction_info";
        open (my $ofh, ">$junction_coords_info_file") or die "Error, cannot write to file: $junction_coords_info_file";
        
        my @fusion_j = reverse sort { $#{$fusion_junctions{$a}} <=> $#{$fusion_junctions{$b}} } keys %fusion_junctions;
        foreach my $fusion (@fusion_j) {
            my @reads = @{$fusion_junctions{$fusion}};
            my $num_fusions = scalar @reads;

            print $ofh "$fusion\t$num_fusions\t" . join(",", @reads) . "\n";
        }
        close $ofh;
    }
    

    exit(0);
}



####
sub get_junction_coord_token {
    my ($scaffold, $original_exon_coord_mapping_href, $genes_matched_aref) = @_;

    my @genes_matched = sort { $a->{coord} <=> $b->{coord} } @$genes_matched_aref;

    my $left_gene = $genes_matched[0]->{gene_id};
    my $right_gene = $genes_matched[$#genes_matched]->{gene_id};
    if ($left_gene eq $right_gene) {
        confess "Error, single gene mapped instead of two";
    }

    my @left_gene_coords;
    my @right_gene_coords;

    foreach my $gene_struct (@genes_matched) {
        my $coord = $gene_struct->{coord};
        my $gene_id = $gene_struct->{gene_id};
        if ($gene_id eq $left_gene) {
            push (@left_gene_coords, $coord);
        }
        else {
            push (@right_gene_coords, $coord);
        }
    }

    @left_gene_coords = sort {$a<=>$b} @left_gene_coords;
    @right_gene_coords = sort {$a<=>$b} @right_gene_coords;

    my $left_fusion_coord = pop @left_gene_coords;
    my $right_fusion_coord = shift @right_gene_coords;

    my $orig_left_fusion_coord = $original_exon_coord_mapping_href->{"$scaffold:$left_fusion_coord"} or die "Error, no coordinate mapping for $scaffold:$left_fusion_coord";
    my $orig_right_fusion_coord = $original_exon_coord_mapping_href->{"$scaffold:$right_fusion_coord"} or die "Error, no coordinate mapping for $scaffold:$right_fusion_coord";

    my $fusion_token = join("\t", $left_gene, $left_fusion_coord, $orig_left_fusion_coord, 
                            $right_gene, $right_fusion_coord, $orig_right_fusion_coord);
    
    return($fusion_token);
}




####
sub hits_exon_bound {
    my ($genome_coords_aref, $exon_bounds_href) = @_;

    my %gene_bounds_matched;
    
    foreach my $coordset (@$genome_coords_aref) {
        my ($lend, $rend) = @$coordset;
        if (my $gene_id = $exon_bounds_href->{$lend}) {
            
            
            $gene_bounds_matched{$gene_id} = { coord => $lend,
                                               gene_id => $gene_id };
        }
        elsif ($gene_id = $exon_bounds_href->{$rend}) {
            
            
            $gene_bounds_matched{$gene_id} = { coord => $rend,
                                               gene_id => $gene_id };
        }

    }
    
    my @genes_matched = values %gene_bounds_matched;

    return(@genes_matched);
    
}


####
sub parse_gtf_file {
    my ($gtf_file, $exon_bounds_href, $original_exon_coord_mapping_href) = @_;

    my %scaff_to_gene_to_coords;

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $scaffold_id = $x[0];
        my $type = $x[2];
        
        unless ($type eq 'exon') { next; }
        
        my $info = $x[8];
        $info =~ /gene_name \"([^\"]+)\"/ or die "Error, cannot parse gene_name from $info";
        my $gene_id = $1;

        my ($lend, $rend) = ($x[3], $x[4]);
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, $lend, $rend);
        $exon_bounds_href->{$scaffold_id}->{$lend} = $gene_id;
        $exon_bounds_href->{$scaffold_id}->{$rend} = $gene_id;
        
        # get original coordinate mapping info
        $info =~ /orig_coord_info \"(\w+),(\d+),(\d+),([+-])\"/ or die "Error, cannot parse original coordinate info from $info";
        my $orig_chr = $1;
        my $orig_lend = $2;
        my $orig_rend = $3;
        my $orig_orient = $4;
        
        my ($orig_end5, $orig_end3) = ($orig_orient eq '+') ? ($orig_lend, $orig_rend) : ($orig_rend, $orig_lend);
        
        $original_exon_coord_mapping_href->{"$scaffold_id:$lend"} = "$orig_chr:$orig_end5";
        $original_exon_coord_mapping_href->{"$scaffold_id:$rend"} = "$orig_chr:$orig_end3";
        
    }
    close $fh;

    
    my %scaffold_to_gene_structs;

    foreach my $scaffold (keys %scaff_to_gene_to_coords) {
        my @genes = keys %{$scaff_to_gene_to_coords{$scaffold}};
    
        my @gene_structs;
    
        foreach my $gene (@genes) {
            my @coords = sort {$a<=>$b} @{$scaff_to_gene_to_coords{$scaffold}->{$gene}};
            my $lend = shift @coords;
            my $rend = pop @coords;
            push (@{$scaffold_to_gene_structs{$scaffold}}, { gene_id => $gene,
                                                             lend => $lend,
                                                             rend => $rend,
                  });
        }
        
    }
        
    return(%scaffold_to_gene_structs);
}


