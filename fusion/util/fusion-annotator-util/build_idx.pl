#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use lib ($ENV{EUK_MODULES});
use TiedHash;


main: {


    my %fusion_to_annots;
    
    { 
        # blast hits
        print STDERR "-parsing blast hits\n";
        open (my $fh, "gunzip -c blastn.gene_pairs.gz | ") or die $!;
        while (<$fh>) {
            chomp;
            my ($geneA, $geneB) = split(/\s+/);
            if ($geneA eq $geneB) { next; }
            $fusion_to_annots{"$geneA--$geneB"}->{"BLASTMATCH"} = 1;
            $fusion_to_annots{"$geneB--$geneA"}->{"BLASTMATCH"} = 1;
        }
        close $fh;
    }


    {
        print STDERR "-parsing prot clusters\n";

        # protein clusters:
        open (my $fh, "gencode.prot_clusters") or die "Error, cannot open file gencode.prot_clusters";
        while (<$fh>) {
            chomp;
            my @genes = split(/\s+/);
            for (my $i = 0; $i < $#genes; $i++) {
                my $gene_i = $genes[$i];
                
                for (my $j = $i + 1; $j <= $#genes; $j++) {
                    my $gene_j = $genes[$j];

                    $fusion_to_annots{"$gene_i--$gene_j"}->{"PROTCLUST"} = 1;
                    $fusion_to_annots{"$gene_j--$gene_i"}->{"PROTCLUST"} = 1;
                    
                }
            }
        }
    }

    {
        
        print STDERR "-parsing nucl clusters\n";
        
        # nucleotide clusters:
        open (my $fh, "gencode.nucl_clusters") or die "Error, cannot open file gencode.nucl_clusters";
        while (<$fh>) {
            chomp;
            my @genes = split(/\s+/);
            for (my $i = 0; $i < $#genes; $i++) {
                my $gene_i = $genes[$i];
                
                for (my $j = $i + 1; $j <= $#genes; $j++) {
                    my $gene_j = $genes[$j];

                    $fusion_to_annots{"$gene_i--$gene_j"}->{"NUCLCLUST"} = 1;
                    $fusion_to_annots{"$gene_j--$gene_i"}->{"NUCLCLUST"} = 1;
                    
                }
            }
        }
        close $fh;
    }
    

    {

        print STDERR "-Parsing fusions found among GTEx normals\n";

        #open (my $fh, "/broad/sigmascratch/WU_CLL_PROJECT/__other_published_studies/Stransky/__SUMMARY/stransky_normals.fusion_to_sample_count") or die $!;
        
        my $stransky_gtex_normals = "Stransky2014.GTEx_normals.pairs";
        open (my $fh, $stransky_gtex_normals) or confess "Error, cannot open file $stransky_gtex_normals";
        
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my ($fusion, $num_normals, $pct_normals) = split(/\t/);
            
            $fusion_to_annots{$fusion}->{"GTEx:$pct_normals%"} = 1;
        }
        close $fh;
    }

    
    {

        print STDERR "- parsing TCGA recurrent fusions from Yoshihara Oncogene 2014\n";
        
        my $yoshihara_tcga_file = "Yoshihara.Oncogene2014.TCGA_recurrent_fusions";
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
            
            my @fusion_annots;
            for (my $i = 0; $i <= 12; $i++) {
                my $num_samples = $x[$i];
                if ($num_samples) {
                    my $tumor = $tumor_type[$i];
                    push (@fusion_annots, "$tumor:$num_samples");
                }
            }
            my $fusion_annot = "TCGA_num_samples[" . join("|", @fusion_annots) . "]";
            
            $fusion_to_annots{$fusion}->{$fusion_annot} = 1;
            
        }
        close $fh;

    }


    {
        # gene spans

        print STDERR "- parsing gene chr locations\n";

        my $gene_spans_file = "gene_spans.txt";
        open (my $fh, $gene_spans_file) or confess "Error, cannot open file $gene_spans_file";
        
        while (<$fh>) {
            chomp;
            
            my ($gene_id, $chr, $lend, $rend, $orient, $gene_name) = split(/\t/);
            if ($gene_name && $gene_name ne ".") {
                $gene_id = $gene_name;
            }
            
            $fusion_to_annots{$gene_id}->{"$chr:$lend-$rend"} = 1;
        }
        close $fh;
        
    }




    {
        # create the database
        
        print STDERR "-writing db idx\n\n";
        
        my $idx = new TiedHash( { create => "fusion_annots.idx" } );

        foreach my $gene_pair (keys %fusion_to_annots) {
            
            my @annots = keys %{$fusion_to_annots{$gene_pair}};
            
            @annots = sort @annots;
            my $annot_string = join(",", @annots);

            $idx->store_key_value($gene_pair, $annot_string);

        }
    }

    print STDERR "-done.\n\n";
    

    exit(0);




}

