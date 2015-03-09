#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 genePairContig.gtf trinity_gmap.gff3\n\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $gff3_align_file = $ARGV[1] or die $usage;

main: {
    

    my %scaffold_to_gene_coordsets = &parse_gtf_file($gtf_file);
    
    #print STDERR Dumper(\%scaffold_to_gene_coordsets);
    
    my %scaffold_to_trinity_coords = &parse_gff3_file($gff3_align_file);
    
    my %trinity_fusion_trans_ids;

    foreach my $scaffold (keys %scaffold_to_gene_coordsets) {
        
        my @genes = keys %{$scaffold_to_gene_coordsets{$scaffold}};

        if (scalar @genes != 2) {
            die "Error, dont have only two genes for scaffold: $scaffold: " . Dumper(\@genes);
        }

        my ($geneA_coords_href, $geneB_coords_href) = &get_gene_coords($scaffold, $scaffold_to_gene_coordsets{$scaffold});
 
        #print "GeneA: " . Dumper($geneA_coords_href) 
        #    . "GeneB: " . Dumper($geneB_coords_href);
        
       
        my @trin_accs = keys (%{$scaffold_to_trinity_coords{$scaffold}});
        foreach my $trin_acc (@trin_accs) {
            my $trin_coords_href = $scaffold_to_trinity_coords{$scaffold}->{$trin_acc};
            
            # ignore singletons
            if (scalar (keys %$trin_coords_href) < 4) { next; } # at least 2 sets of coordinates, indicating an intron
            
            #print "Trinity: $trin_acc " . Dumper($trin_coords_href);

            if (&shared_coordinate($geneA_coords_href, $trin_coords_href)
                &&
                &shared_coordinate($geneB_coords_href, $trin_coords_href) ) {

                my ($break_left, $break_right) = &get_breakpoint_coords($geneA_coords_href, $geneB_coords_href, $trin_coords_href);

                $trinity_fusion_trans_ids{$trin_acc} = "$scaffold:$break_left-$break_right";
            }

        }
    }

    &report_trin_fusions($gff3_align_file, \%trinity_fusion_trans_ids);
    
    

    exit(0);
}


####
sub shared_coordinate {
    my ($coordsA_href, $coordsB_href) = @_;

    foreach my $coord (keys %$coordsA_href) {
        
        if ($coordsB_href->{$coord}) {
            return(1);
        }
    }

    return(0);
}



####
sub get_gene_coords {
    my ($scaffold, $genes_to_coords_href) = @_;

    my @gene_coords_hrefs;
    foreach my $gene (split(/--/, $scaffold)) {
        
        my @coords = @{$genes_to_coords_href->{$gene}};

        my %gene_coords;
        foreach my $coordpair (@coords) {
            my ($lend, $rend) = @$coordpair;
            $gene_coords{$lend} = 1;
            $gene_coords{$rend} = 1;
        }
        push (@gene_coords_hrefs, \%gene_coords);
    }

    return(@gene_coords_hrefs);

}


####
sub parse_gtf_file {
    my ($gtf_file) = @_;

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
        push (@{$scaff_to_gene_to_coords{$scaffold_id}->{$gene_id}}, [$lend, $rend]);
        
    }
    close $fh;

    
    return(%scaff_to_gene_to_coords);
}




####
sub parse_gff3_file {
    my ($gff3_align_file) = @_;

    
    my %scaffold_to_trans_coords;
    
    open (my $fh, $gff3_align_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; } # comment line
        unless (/\w/) { next; }
        
        chomp;
        my @x = split(/\t/);
        my $scaff = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $info = $x[8];

        my $trinity_id;
        if ($info =~ /ID=([^;]+)\.path[\d+];/) {

            $trinity_id = $1;
        }
        else {
            die "Error, cannot find Trinity ID from $info";
        }

        $scaffold_to_trans_coords{$scaff}->{$trinity_id}->{$lend} = 1;
        $scaffold_to_trans_coords{$scaff}->{$trinity_id}->{$rend} = 1;
                
    }
    close $fh;
    
    return(%scaffold_to_trans_coords);

}

####
sub report_trin_fusions {
    my ($gff3_align_file, $trin_ids_href) = @_;

    
    foreach my $trin_id (keys %$trin_ids_href) {
        my $scaff_breakpoint = $trin_ids_href->{$trin_id};
        print "#TrinityFusionTranscript:\t$trin_id\t$scaff_breakpoint\n";
    }
        
    my %scaffold_to_trans_coords;
    
    open (my $fh, $gff3_align_file) or die $!;
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        my $line = $_;
        my @x = split(/\t/);
        my $info = $x[8];


        if ($info =~ /ID=([^;]+)\.path[\d+];/) {

            my $trinity_id = $1;
            
            if ($trin_ids_href->{$trinity_id}) {
                print $line;
                
            }
        }
        else {
            die "Error, cannot find Trinity ID from $info";
        }
        
    }
    close $fh;


    return;
}

####
sub get_breakpoint_coords {
    my ($geneA_coords_href, $geneB_coords_href, $trin_coords_href) = @_;


    ## get left breakpoint
    my @left_shared_coords;

    foreach my $coord (keys %$geneA_coords_href) {
        if ($trin_coords_href->{$coord}) {
            push (@left_shared_coords, $coord);
        }
    }
    
    @left_shared_coords = sort {$a<=>$b} @left_shared_coords;
    my $left_breakpoint = pop @left_shared_coords;
    unless ($left_breakpoint) {
        confess "Error, no left breakpoint";
    }
    
    ## get right breakpoint
    my @right_shared_coords;

    foreach my $coord (keys %$geneB_coords_href) {
        if ($trin_coords_href->{$coord}) {
            push (@right_shared_coords, $coord);
        }
    }
    @right_shared_coords = sort {$a<=>$b} @right_shared_coords;

    my $right_breakpoint = shift @right_shared_coords;

    unless ($right_breakpoint) {
        confess "Error, no right breakpoint";
    }
       
    return($left_breakpoint, $right_breakpoint);
}

