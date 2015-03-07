#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use __GLOBALS__;
use Set::IntervalTree;
use SAM_reader;
use SAM_entry;
use FusionAnnotator;


my $usage = "usage: $0 gsnap.split.sam, ...\n\n";

my @sam_files = @ARGV;
unless (@sam_files) {
    die $usage;
}

my $gtf_file = "/seq/regev_genome_portal/RESOURCES/human/Hg19/gencode.v19/gencode.v19.rna_seq_pipeline.gtf.exons";


main: {
 
    ## prep hg19 resource data
    my %chr_to_interval_tree;
    my %chr_coord_to_splice_junction;


    print STDERR "-parsing annotation coordinate info.\n";
    &parse_junctions($gtf_file, \%chr_to_interval_tree, \%chr_coord_to_splice_junction);


    my @to_genes_files;
    
    foreach my $sam_file (@sam_files) {

        print STDERR "-examining file: $sam_file\n";

        my $sam_reader = new SAM_reader($sam_file);

        my $to_genes_file = "$sam_file.to_genes";
        open (my $ofh, ">$to_genes_file") or die "Error, cannot write to file $to_genes_file";
        
        push (@to_genes_files, $to_genes_file);
        
        while (my $sam_entry = $sam_reader->get_next()) {
            
            my $core_read_name = $sam_entry->get_core_read_name();
            
            my $chr = $sam_entry->get_scaffold_name();
            
            unless (exists $chr_to_interval_tree{$chr}) { next; } # no annotations to search overlaps of.
            
            my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
            
            my @coords;
            my %overlapping_genes;
            foreach my $exon_segment (@$genome_coords_aref) {
                my ($lend, $rend) = @$exon_segment;
                push (@coords, $lend, $rend);
                my $overlaps_aref = $chr_to_interval_tree{$chr}->fetch($lend, $rend);
                if (@$overlaps_aref) {
                    foreach my $overlapping_gene (@$overlaps_aref) {
                        $overlapping_genes{$overlapping_gene} = 1;
                    }
                }
            }
            
            # check for splice junction 
            @coords = sort {$a<=>$b} @coords;

            my $lend = shift @coords;
            my $rend = shift @coords;
            
            if (my $splice_bound = $chr_coord_to_splice_junction{"$chr:$lend"}) {
                $overlapping_genes{$splice_bound} = 1;
            }
            if (my $splice_bound = $chr_coord_to_splice_junction{"$chr:$rend"}) {
                $overlapping_genes{$splice_bound} = 1;
            }

            if (%overlapping_genes) {
                print $ofh "$core_read_name\t" . join(",", sort keys %overlapping_genes) . "\n";
            }
                        
        } # end for sam entry

        close $ofh;

    } # end for sam file

    print STDERR "-identifying candidate fusions.\n";

    my $sorted_candidates_file = "gsnap.sorted_candidates.tmp";
    my $cmd = "sort -k1,1 @to_genes_files > $sorted_candidates_file";
    &process_cmd($cmd);

    report:

    $sorted_candidates_file = "gsnap.sorted_candidates.tmp";
    print STDERR "-retrieving list of fusion predictions.\n";
    my $gsnap_fusion_candidates = "gsnap.fusions.list\n";
    &identify_gsnap_fusions($sorted_candidates_file, $gsnap_fusion_candidates);
    
    
    exit(0);
    
}


####
sub identify_gsnap_fusions {
    my ($candidates_input_file, $fusion_list_output_file) = @_;
    
    my $prev_read_name = "";
    my %indiv_fragment_fusion_genes;

    my %fusion_candidates;

    open (my $fh, $candidates_input_file) or die "Error, cannot open file $candidates_input_file";
    while (<$fh>) {
        chomp;
        my ($read_name, $gene_list) = split(/\t/);
        if ($read_name ne $prev_read_name) {
            if (%indiv_fragment_fusion_genes) {
                &examine_fusion_candidates(\%fusion_candidates, \%indiv_fragment_fusion_genes);
            }
            %indiv_fragment_fusion_genes = (); #re-init
        }
        
        foreach my $gene (split(/,/, $gene_list) ) {
            $indiv_fragment_fusion_genes{$gene}++;
        }
        
        $prev_read_name = $read_name;
    }
    
    # get last one
    if (%indiv_fragment_fusion_genes) {
        &examine_fusion_candidates(\%fusion_candidates, \%indiv_fragment_fusion_genes);
    }
 
    
  report_fusions: {

      ## output the fusion list, require junction reads for candidates.
      open (my $ofh, ">$fusion_list_output_file") or die "Error, cannot write to file $fusion_list_output_file";
      print $ofh "#fusion_genes\tjunction_reads\tspanning_frags\n";
      if (my $fusion_genes_href = $fusion_candidates{JUNCTION}) {
          my @fusion_genes = keys %{$fusion_genes_href};
          foreach my $fusion_gene (@fusion_genes) {
              my $junction_count = $fusion_genes_href->{$fusion_gene};
              
              my $spanning_gene_key = join("--", sort split(/--/, $fusion_gene));
              my $spanning_count = $fusion_candidates{SPANNING}->{$spanning_gene_key} || 0;
              
              my ($geneA, $geneB) = split(/--/, $fusion_gene);
              
              ## no self-fusions ;)
              if ($geneA eq $geneB) { next; }

              my @annotations = &FusionAnnotator::get_annotations($geneA, $geneB);
              
              print $ofh join("\t", $fusion_gene, $junction_count, $spanning_count, 
                  join(",", @annotations)) . "\n";
          }
      }
      close $ofh;
    }


    print STDERR "-done.  See $fusion_list_output_file for fusion candidates.\n\n";

    exit(0);
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


####
sub examine_fusion_candidates {
    my ($fusion_candidates_href, $indiv_fragment_fusion_genes_href) = @_;

    my %donors;
    my %acceptors;

    my %genes;
    foreach my $gene_info (keys %$indiv_fragment_fusion_genes_href) {

        my ($gene_id, $donor_or_acceptor) = split(/::/, $gene_info);
        
        $genes{$gene_id}++;
        if ($donor_or_acceptor) {
            if ($donor_or_acceptor eq 'DONOR') {
                $donors{$gene_id}++;
            }
            elsif ($donor_or_acceptor eq 'ACCEPTOR') {
                $acceptors{$gene_id}++;
            }
            else {
                die "Error, got donor_or_acceptor: $donor_or_acceptor, but not DONOR or ACCEPTOR";
            }
        }
    }

    my @genes_involved = keys %genes;
    
    my @donor_genes = keys %donors;
    my @acceptor_genes = keys %acceptors;
    
    if (scalar @donor_genes == 1 && scalar @acceptor_genes == 1) {
        ## got directional fusion w/ junction read
        my $fusion_name = $donor_genes[0] . "--" . $acceptor_genes[0];
        $fusion_candidates_href->{JUNCTION}->{$fusion_name}++;
    }
    elsif (scalar @genes_involved == 2) {
        ## consider it a spanning read
        @genes_involved = sort @genes_involved;
        my $fusion_name = join("--", @genes_involved);
        $fusion_candidates_href->{SPANNING}->{$fusion_name}++;
    }
    
    return;
}
            



####
sub describe_read_pair_info {
    my ($read_name, $read_info_aref) = @_;

    my @read_1_genes;
    my @read_2_genes;

    my $read_1_genes_href = $read_info_aref->[1];
    if ($read_1_genes_href) {
        @read_1_genes = sort keys (%$read_1_genes_href);
    }

    my $read_2_genes_href = $read_info_aref->[2];
    if ($read_2_genes_href) {
        @read_2_genes = sort keys (%$read_2_genes_href);
    }
    
    if (@read_1_genes && @read_2_genes) {
        my $fusion_token = join(",", @read_1_genes) . "--" . join(",", @read_2_genes);
                
        print "$read_name\t$fusion_token\n";
    }
    
    return;
}


####
sub add_overlapping_genes {
    my ($read_info_aref, $overlaps_aref, $read_dir) = @_;
    
    foreach my $gene (@$overlaps_aref) {
        $read_info_aref->[$read_dir]->{$gene} = 1;
    }

    return;
}


####
sub parse_junctions {
    my ($gtf_file, $chr_to_interval_tree_href, $chr_coord_to_splice_junction_href) = @_;


    my $counter = 0;
    
    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
                
        my @x = split(/\t/);
        my $chr = $x[0];
        my $type = $x[2];
    
        unless ($type eq 'exon') { next; }
        
        $counter++;
        #print STDERR "\r[$counter]    " if $counter % 1000 == 0;
        

        my $lend = $x[3];
        my $rend = $x[4];
        
        my $orient = $x[6];


                
        my $info = $x[8];

        my $gene_id;
        if ($info =~ /gene_id \"?([^;\"]+)\"?;/) {
            $gene_id = $1;
        }

        if ($info =~ /gene_name \"?([^;\"]+)\"?;/) {
            $gene_id = $1;
        }

        unless ($gene_id) {
            confess "Error, cannot find gene_id from $info";
        }

        if ($rend - $lend > 1) {
            &add_chr_interval($chr_to_interval_tree_href, $chr, $gene_id, $lend, $rend);
            
            &add_splice_junction($chr_coord_to_splice_junction_href, $chr, $gene_id, $lend, $rend, $orient);
        }

    }
    
    return;

}
        

####
sub add_splice_junction {
    my ($chr_coord_to_splice_junction_href, $chr, $gene_id, $lend, $rend, $orient) = @_;

    my ($donor, $acceptor) = ($orient eq '+') ? ($rend, $lend) : ($lend, $rend);

    $chr_coord_to_splice_junction_href->{"$chr:$donor"} = join("::", $gene_id, "DONOR");
    $chr_coord_to_splice_junction_href->{"$chr:$acceptor"} = join("::", $gene_id, "ACCEPTOR");
    
    return;
}




####
sub add_chr_interval {
    my ($chr_exons_interval_tree_href, $chr, $gene_id, $lend, $rend) = @_;

    my $i_tree = $chr_exons_interval_tree_href->{$chr};
    unless ($i_tree) {
        $i_tree = $chr_exons_interval_tree_href->{$chr} = Set::IntervalTree->new;
    }
    
    eval {
        $i_tree->insert($gene_id, $lend, $rend);
    };

    if ($@) {
        print STDERR "Error, $@, $gene_id,$lend,$rend";
    }
    return;
}
        
