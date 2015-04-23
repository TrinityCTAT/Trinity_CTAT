#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 genePairContig.gtf read_alignments.bam bam.fusion_junction_info\n\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $bam_file = $ARGV[1] or die $usage;
my $junction_info_file = $ARGV[2] or die $usage; # ignore these among the spanning reads

main: {
    
    my %junction_reads_ignore;
    my %fusion_junctions;
    my %fusion_breakpoint_info;
    
    {
        open (my $fh, "$junction_info_file") or die "error, cannot open file $junction_info_file";
        while (<$fh>) {
            chomp;
            my ($geneA, $coordA, $orig_coordA, $geneB, $coordB, $orig_coordB, $splice_info, $fusion_read_count, $fusion_read_list) = split(/\t/);

            foreach my $fusion_read (split(/,/, $fusion_read_list)) {
                $junction_reads_ignore{$fusion_read}++;
            }
            
            my $fusion_name = join("--", $geneA, $geneB);
            push (@{$fusion_junctions{$fusion_name}}, "$coordA-$coordB");
            
            my $breakpoint = "$coordA-$coordB";
            
            $fusion_breakpoint_info{"$fusion_name|$breakpoint"} = join("\t", $geneA, $coordA, $orig_coordA,
                                                                       $geneB, $coordB, $orig_coordB, $splice_info);
            

        }
        close $fh;
    }
    

    my %exon_bounds;
    my %scaffold_to_gene_structs = &parse_gtf_file($gtf_file, \%exon_bounds);
    
    my %scaffold_to_gene_breaks;
    {
        foreach my $scaffold (keys %scaffold_to_gene_structs) {
            my @gene_structs = @{$scaffold_to_gene_structs{$scaffold}};
            if (scalar @gene_structs != 2) {
                die "Error, didn't find only 2 genes in the gtf file: " . Dumper(\@gene_structs);
            }
            
            @gene_structs = sort {$a->{lend} <=> $b->{lend}} @gene_structs;
            
            my $left_gene = $gene_structs[0];
            my $right_gene = $gene_structs[1];
            
            my $gene_bound_left = $left_gene->{rend};
            my $gene_bound_right = $right_gene->{lend};
            
            
            if ($gene_bound_left > $gene_bound_right) {
                die "Error, bounds out of order: $gene_bound_left ! <  $gene_bound_right";
            }
            $scaffold_to_gene_breaks{$scaffold} = [$gene_bound_left, $gene_bound_right];
        }
    }
    
    print STDERR "Scaffold to gene breaks: " . Dumper(\%scaffold_to_gene_breaks);
    

    ## for each paired read, get the bounds of that read
    my %scaffold_read_pair_to_read_bounds;    
    my %core_counter;
    {
        ## find the reads that matter:
        my $counter = 0;
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            $counter++;
            print STDERR "\r[$counter]   " if $counter % 1000 == 0;
            
            #print STDERR Dumper($sam_entry);

            my $qual_val = $sam_entry->get_mapping_quality();
            #unless ($qual_val >= $MIN_QUALITY) { next; }
            

            my $scaffold = $sam_entry->get_scaffold_name();
            my $read_name = $sam_entry->get_read_name();
            my $token = join("$;", $read_name, $scaffold);

            my $full_read_name = $sam_entry->reconstruct_full_read_name();            

            if ($junction_reads_ignore{$full_read_name}) { next; }
            
        
            my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
            unless (&overlaps_exon($genome_coords_aref, $exon_bounds{$scaffold}) ) { 
                #print STDERR "No exon overlap: " . Dumper($genome_coords_aref) . Dumper($exon_bounds{$scaffold});
                next; 
            } # only examine exon-overlapping entries
            
            #print STDERR "got exon overlap: " . Dumper($genome_coords_aref) . Dumper($exon_bounds{$scaffold});
            
            my ($span_lend, $span_rend) = sort {$a<=>$b} $sam_entry->get_genome_span();

            
            
            
            #print STDERR "$span_lend to $span_rend\n";
            
            if ($full_read_name =~ /^(\S+)\/([12])$/) {
                my ($core, $pair_end) = ($1, $2);
                $core_counter{"$scaffold|$core"}++;  # track how many alignments we have for this rnaseq fragment
                
                $scaffold_read_pair_to_read_bounds{$scaffold}->{$core}->[$pair_end-1] = [$span_lend, $span_rend];
                
            }
        }
    }
   
    #print STDERR "Read pair to read bounds: " . Dumper(\%scaffold_read_pair_to_read_bounds);
    
    
    my %fusion_to_spanning_reads;

    my %fusion_to_contrary_support;
    
    # determine which reads are spanning reads
    my %spanning_read_want;
    {
        foreach my $scaffold (keys %scaffold_to_gene_breaks) {
            my ($gene_bound_left, $gene_bound_right) = @{$scaffold_to_gene_breaks{$scaffold}};
            
            if ($gene_bound_left > $gene_bound_right) { 
                die "Error, gene bounds out of range for $scaffold: $gene_bound_left - $gene_bound_right "; 
            }
            
            foreach my $fragment (keys %{$scaffold_read_pair_to_read_bounds{$scaffold}}) {

                if ($core_counter{"$scaffold|$fragment"} != 2) { next; } # ignore those fragments that have multiply-mapping reads to this contig.
                my @pair_coords = grep { defined $_ } @{$scaffold_read_pair_to_read_bounds{$scaffold}->{$fragment}};
                if (scalar @pair_coords > 1) {
                    # need both paired ends
                    
                    @pair_coords = sort {$a->[0] <=> $b->[0]} @pair_coords;
                    
                    my $left_read_rend = $pair_coords[0]->[1];
                    my $right_read_lend = $pair_coords[1]->[0];

                    foreach my $fusion (@{$fusion_junctions{$scaffold}}) {
                        my ($break_lend, $break_rend) = split(/-/, $fusion);
                        
                        if ($left_read_rend < $break_lend && $right_read_lend > $break_rend) {
                            
                            $spanning_read_want{"$scaffold|$fragment"}++; # capture for SAM-retreival next.
                        
                            push (@{$fusion_to_spanning_reads{"$scaffold|$fusion"}}, $fragment);
                        }
                        elsif ($left_read_rend < $break_lend && $break_lend < $right_read_lend
                               && $right_read_lend < $gene_bound_right) {

                            ## contrary support at left junction
                            push (@{$fusion_to_contrary_support{"$scaffold|$fusion"}->{left}}, $fragment);
                        }
                        elsif ($left_read_rend < $break_rend && $break_rend < $right_read_lend
                               && $left_read_rend > $gene_bound_left) {
                            ## contrary support at right junction
                            push (@{$fusion_to_contrary_support{"$scaffold|$fusion"}->{right}}, $fragment);
                        }
                        
                    }
                }
            }
        }
    }
    
    #print STDERR "spanning reads want: " . Dumper(\%spanning_read_want);
    

    # output the spanning reads we want
    if (%spanning_read_want) {
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            
            my $qual_val = $sam_entry->get_mapping_quality();
            #unless ($qual_val >= $MIN_QUALITY) { next; }
            
            
            my $scaffold = $sam_entry->get_scaffold_name();
            my $core_read_name = $sam_entry->get_core_read_name();
            if ($spanning_read_want{"$scaffold|$core_read_name"}) {
                print $sam_entry->get_original_line() . "\n";
            }
        }
    }
    
    {
        
        print STDERR "-outputting the spanning read info.\n";

        ## output the spanning read info
        my $spanning_read_info_file = "$bam_file.fusion_spanning_info";
        open (my $ofh, ">$spanning_read_info_file") or die "Error, cannot write to $spanning_read_info_file";
        foreach my $fusion_n_breakpoint (sort keys %fusion_to_spanning_reads) {
            my ($fusion_name, $breakpoint) = split(/\|/, $fusion_n_breakpoint);
            my ($geneA, $geneB) = split(/--/, $fusion_name);

            my $fusion_info = $fusion_breakpoint_info{$fusion_n_breakpoint};
            
            my ($coordA, $coordB) = split(/-/, $breakpoint);
            
            my @spanning_reads = @{$fusion_to_spanning_reads{$fusion_n_breakpoint}};
            
            my $num_spanning = scalar(@spanning_reads);
            
            ## examine contrary support
            my @contrary_left_support;
            if (my $support_aref = $fusion_to_contrary_support{$fusion_n_breakpoint}->{left}) {
                @contrary_left_support = @$support_aref;
            }
            my @contrary_right_support;
            if (my $support_aref = $fusion_to_contrary_support{$fusion_n_breakpoint}->{right}) {
                @contrary_right_support = @$support_aref;
            }
            
            
            my $num_left_contrary_support = scalar(@contrary_left_support);
            my $num_right_contrary_support = scalar(@contrary_right_support);
            

            print $ofh join("\t", $fusion_info, $num_spanning) . "\t" . join(",", @spanning_reads)
                . "\t$num_left_contrary_support\t" . join(",", @contrary_left_support) 
                . "\t$num_right_contrary_support\t" . join(",", @contrary_right_support)
                . "\n";
        }
        close $ofh;
    }
    
    
    exit(0);
}

####
sub overlaps_exon {
    my ($genome_coords_aref, $exon_bounds_href) = @_;

    foreach my $coordset (@$genome_coords_aref) {
        my ($lend, $rend) = @$coordset;
        foreach my $exon_coordset (keys %$exon_bounds_href) {
            my ($e_lend, $e_rend) = split(/-/, $exon_coordset);
            if ($e_lend < $rend && $e_rend > $lend) {
                return(1);
            }
        }
    }
    return(0); # no dice
}



####
sub parse_gtf_file {
    my ($gtf_file, $exon_bounds_href) = @_;

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
        $exon_bounds_href->{$scaffold_id}->{"$lend-$rend"} = 1;
        
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


