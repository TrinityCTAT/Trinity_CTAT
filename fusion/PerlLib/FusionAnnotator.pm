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

    if (my $dist_annot = &__get_distance_annotation($geneA, $geneB)) {
        push (@annotations, $dist_annot);
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
        return(undef);
    }


    #print STDERR "A: $chr_info_A\tB: $chr_info_B\n";
    

    my ($chrA, $coords_A) = split(/:/, $chr_info_A);
    $coords_A =~ s/\,.*$//;
    my ($lendA, $rendA) = split(/-/, $coords_A);
    
    my ($chrB, $coords_B) = split(/:/, $chr_info_B);
    $coords_B =~ s/\,.*$//;
    my ($lendB, $rendB) = split(/-/, $coords_B);
    
    if ($chrA ne $chrB) {
        return(undef);
    }

    my $dist = -1;
    if ($lendA < $rendB && $rendA > $lendB) {
        # overlap
        $dist = 0;
    }
    my @coords = sort {$a<=>$b} ($lendA, $rendA, $lendB, $rendB);
    $dist = $coords[2] - $coords[1];

    if ($dist > 0 && $dist <= $MAX_NEIGHBOR_DIST) {
        return("NEIGHBORS\[$dist]");
    }
    else {
        return(undef);
    }
}


    




sub init {
    
    $IDX = new TiedHash( { use => "$FUSION_ANNOTATOR_LIB/fusion_annots.idx" } );
    
}

1; #EOM


    
        
                
      
