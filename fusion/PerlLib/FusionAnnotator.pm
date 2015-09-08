package FusionAnnotator;

use strict;
use warnings;
use Carp;

BEGIN {
    unless ($ENV{FUSION_ANNOTATOR_LIB}) {
        confess "ERROR, must have environmental variable FUSION_ANNOTATOR_LIB set to the path to find the fusion annotator module and data files.";
    }
}

use FindBin;
use lib ("$FindBin::Bin");
use TiedHash;


my $MAX_NEIGHBOR_DIST = 100000;

my $FUSION_ANNOTATOR_LIB = $ENV{FUSION_ANNOTATOR_LIB};


my $IDX; # tied hash

####
sub set_max_neighbor_dist {
    my ($dist) = @_;

    $MAX_NEIGHBOR_DIST = $dist;
    
    return;
}

####
sub get_annotations {
    my ($geneA, $geneB) = @_;
    
    unless ($IDX) {
        &init();
    }
    
    my $annotation_text = $IDX->get_value("$geneA--$geneB");
    
    my @annotations;
    if ($annotation_text) {
        push (@annotations, $annotation_text);
    }

    if (my @dist_annots = &__get_distance_annotation($geneA, $geneB)) {
        push (@annotations, @dist_annots);
    }
    
    
    return(@annotations);

}


####
sub __get_distance_annotation {
    my ($geneA, $geneB) = @_;

    my $chr_info_A = $IDX->get_value($geneA);
    my $chr_info_B = $IDX->get_value($geneB);
    
    unless ($chr_info_A && $chr_info_B) {
        # cant compare them
        return();
    }


    #print STDERR "A: $chr_info_A\tB: $chr_info_B\n";
    

    my ($chrA, $coords_A, $orientA) = split(/:/, $chr_info_A);
    $coords_A =~ s/\,.*$//;
    my ($lendA, $rendA) = split(/-/, $coords_A);
    
    my ($chrB, $coords_B, $orientB) = split(/:/, $chr_info_B);
    $coords_B =~ s/\,.*$//;
    my ($lendB, $rendB) = split(/-/, $coords_B);
    
    my $dist = -1;
    if ($lendA < $rendB && $rendA > $lendB) {
        # overlap
        $dist = 0;
    }

    my @annotations;

    if ($chrA eq 'chrM' || $chrB eq 'chrM') {
        push (@annotations, "INCL_MITOC_GENE");
    }

    if ($chrA eq $chrB) {
    
        my @coords = sort {$a<=>$b} ($lendA, $rendA, $lendB, $rendB);
        $dist = $coords[2] - $coords[1];
        
        
        
        if ($dist > 0 && $dist <= $MAX_NEIGHBOR_DIST) {
            
            if ($lendA < $rendB && $rendA > $lendB) {
                push (@annotations, "NEIGHBORS_OVERLAP:$orientA:$orientB:[$dist]");
            }
            elsif ($orientA ne $orientB) { 
                push(@annotations, "LOCAL_INVERSION:$orientA:$orientB:[$dist]");
            }
            elsif ($orientA eq '+' && $lendB < $rendA) { 
                push (@annotations, "LOCAL_REARRANGEMENT:$orientA:[$dist]");
            }
            elsif ($orientA eq '-' && $rendB > $lendA) { 
                push (@annotations, "LOCAL_REARRANGEMENT:$orientA:[$dist]"); 
            }
            else {
                # no other weirdness, just neighbors, probably readthru transcription
                
                push (@annotations, "NEIGHBORS\[$dist]");
            }
        }
    }

    
    return(@annotations);
}


    




sub init {
    
    $IDX = new TiedHash( { use => "$FUSION_ANNOTATOR_LIB/fusion_annots.idx" } );
    
}

1; #EOM


    
        
                
      
