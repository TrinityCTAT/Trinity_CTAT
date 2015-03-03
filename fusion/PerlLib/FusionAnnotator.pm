package FusionAnnotator;

use strict;
use warnings;
use Carp;

BEGIN {
    unless ($ENV{FUSION_ANNOTATOR_LIB}) {
        confess "ERROR, must have environmental variable FUSION_ANNOTATOR_LIB set to the path to find the fusion annotator module and data files.";
    }
}



our %blastmatch;
our %paralogs;
our %stransky_GTEx_normals;
our %YoshiharaTCGA;
our %gene_spans;


my $MAX_NEIGHBOR_DIST = 100000;

my $FUSION_ANNOTATOR_LIB = $ENV{FUSION_ANNOTATOR_LIB};

####
sub set_max_neighbor_dist {
    my ($dist) = @_;

    $MAX_NEIGHBOR_DIST = $dist;
    
    return;
}

####
sub get_annotations {
    my ($geneA, $geneB) = @_;
    
    unless (%blastmatch) {
        &init();
    }
    
    my @annotations;

    if ($blastmatch{$geneA}->{$geneB}) {
        push (@annotations, "BLASTMATCH");
    }

    if ($paralogs{$geneA}->{$geneB}) {
        push (@annotations, "PARALOGS");
    }
    
    if (my $stransky_vals = $stransky_GTEx_normals{$geneA}->{$geneB}) {
        push (@annotations, "StranskyGTEx($stransky_vals)");
    }
    
    if (my $yoshihara_vals = $YoshiharaTCGA{$geneA}->{$geneB}) {
        push (@annotations, $yoshihara_vals);
    }
    
    if (my $dist_annot = &__get_distance_annotation($geneA, $geneB)) {
        push (@annotations, $dist_annot);
    }
    

    return(@annotations);
}

####
sub __get_distance_annotation {
    my ($geneA, $geneB) = @_;

    my $chr_info_A = $gene_spans{$geneA};
    my $chr_info_B = $gene_spans{$geneB};
    
    unless ($chr_info_A && $chr_info_B) {
        # cant compare them
        return(undef);
    }

    my ($chrA, $lendA, $rendA) = split(/[:-]/, $chr_info_A);
    my ($chrB, $lendB, $rendB) = split(/[:-]/, $chr_info_B);

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
    
    BlastMatches: {
        print STDERR "FusionAnnotator.pm - parsing gene blast matches.\n";
        
        
        #open (my $fh, "/seq/regev_genome_portal/RESOURCES/human/Hg19/gencode.v19/blast.1e-10.101.outfmt6.wLens.genes") or die $!;
        my $blast_matches_file = "$FUSION_ANNOTATOR_LIB/blast_matches.dat";
        open (my $fh, $blast_matches_file) or confess "Error, cannot open file $blast_matches_file";
        
        while (<$fh>) {
            my @x = split(/\t/);
            my ($gene_i, $gene_j) = ($x[0], $x[1]);
            $blastmatch{$gene_i}->{$gene_j} = 1;
            $blastmatch{$gene_j}->{$gene_i} = 1;
        }
        
        close $fh;
    }
    


    Paralogs: { 

        print STDERR "FusionAnnotator.pm - parsing paralog info.\n";
        
        #open (my $fh, "/seq/regev_genome_portal/RESOURCES/human/Hg19/gencode.v19/pairs.genes.slclust.j0.2.dat") or die $!;
        my $paralog_file = "$FUSION_ANNOTATOR_LIB/paralogs.dat";
        open (my $fh, $paralog_file) or confess "Error, cannot open file $paralog_file";
        
        while (<$fh>) {
            chomp;
            my @x = split(/\s+/);
            for (my $i = 0; $i < $#x; $i++) {
                my $gene_i = $x[$i];
                
                for (my $j = $i + 1; $j <= $#x; $j++) {
                    my $gene_j = $x[$j];
                    
                    $paralogs{$gene_i}->{$gene_j} = 1;
                    
                    $paralogs{$gene_j}->{$gene_i} = 1;
                    
                }
            }
        }
        close $fh;
    }


    Normals: {

        print STDERR "FusionAnnotator.pm - parsing fusions found among GTEx normals\n";

        #open (my $fh, "/broad/sigmascratch/WU_CLL_PROJECT/__other_published_studies/Stransky/__SUMMARY/stransky_normals.fusion_to_sample_count") or die $!;
        
        my $stransky_gtex_normals = "$FUSION_ANNOTATOR_LIB/Stransky2014.GTEx_normals.pairs";
        open (my $fh, $stransky_gtex_normals) or confess "Error, cannot open file $stransky_gtex_normals";
        
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my ($fusion, $num_normals, $pct_normals) = split(/\t/);

            my ($geneA, $geneB) = split(/--/, $fusion);
            $stransky_GTEx_normals{$geneA}->{$geneB} = "$num_normals,$pct_normals\%";
        }
        close $fh;
      }


    Tumors: {

        print STDERR "FusionAnnotator.pm - parsing TCGA recurrent fusions from Yoshihara Oncogene 2014\n";
        
        #open (my $fh, "/broad/sigmascratch/WU_CLL_PROJECT/__other_published_studies/Yoshihara_Oncogene_2014/onc2014406x2.txt") or die $!;
        my $yoshihara_tcga_file = "$FUSION_ANNOTATOR_LIB/Yoshihara.Oncogene2014.TCGA_recurrent_fusions";
        open (my $fh, $yoshihara_tcga_file) or confess "Error, cannot open file $yoshihara_tcga_file";
        
        my $topline = <$fh>;
        my $header = <$fh>;
        chomp $header;
        my @tumor_type = split(/\t/, $header);
        shift @tumor_type;
        
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $fusion = shift @x;
            my ($geneA, $geneB) = split(/__/, $fusion);
            
            my @fusion_annots;
            for (my $i = 0; $i <= 12; $i++) {
                my $num_samples = $x[$i];
                if ($num_samples) {
                    my $tumor = $tumor_type[$i];
                    push (@fusion_annots, "$tumor:$num_samples");
                }
            }
            my $fusion_annot = "YoshiharaTCGA[" . join("|", @fusion_annots) . "]";

            $YoshiharaTCGA{$geneA}->{$geneB} = $fusion_annot;
        }
        close $fh;

      }
                
    gene_spans: {
        
        print STDERR "FusionAnnotator.pm - parsing gene chr locations\n";

        #open (my $fh, "/seq/regev_genome_portal/RESOURCES/human/Hg19/Annotations/Hg19_Gencode/gencode.v19.rna_seq_pipeline.gtf.gene_spans") or die $!;
        my $gene_spans_file = "$FUSION_ANNOTATOR_LIB/gene_spans.txt";
        open (my $fh, $gene_spans_file) or confess "Error, cannot open file $gene_spans_file";
        
        while (<$fh>) {
            chomp;
            
            my ($gene_id, $chr, $lend, $rend, $orient, $gene_name) = split(/\t/);
            if ($gene_name && $gene_name ne ".") {
                $gene_id = $gene_name;
            }

            $gene_spans{$gene_id} = "$chr:$lend-$rend";
        }
        close $fh;
        
      }

      return;

}

1; #EOM


    
        
                
      
