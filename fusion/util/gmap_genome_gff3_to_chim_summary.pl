#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);



my $usage = <<__EOUSAGE__;

#####################################################################
#
# Required:
#
#  --gmap_gff3 <string>        gmap gff3 alignment output
#
####################
#
# Optional:
#
#  --annot_gtf <string>        transcript structures in gtf file format.
#
######################################################################


__EOUSAGE__

    ;


my $help_flag;
my $gmap_gff3_file;
my $annot_gtf_file;

&GetOptions ( 'h' => \$help_flag,
              'gmap_gff3=s' => \$gmap_gff3_file,
              'annot_gtf=s' => \$annot_gtf_file,
              
              );


if ($help_flag) {
    die $usage;
}

unless ($gmap_gff3_file) {
    die $usage;
}


my %genes;

main: {

    if ($annot_gtf_file) {
        open (my $fh, $annot_gtf_file) or die "Error, cannot open file $annot_gtf_file";
        while (<$fh>) {
            chomp;
            unless (/\w/) { next; }
            if (/^\#/) { next; }
            s/^>//;
            my @x = split(/\t/);
            
            unless ($x[2] eq "exon") { next; }
            
            my $info = $x[8];
            $info =~ /gene_id \"([^\"]+)/ or die "Error, cannot extract gene_id from $_ [specifically from: $info]";
            my $gene_id = $1 or die "Error, no gene_id from $_";

            if ($info =~ /gene_name \"([^\"]+)/) {
                # use gene name instead
                $gene_id = $1;
            }
            
            $info =~ /transcript_id \"([^\"]+)/ or die "Error, cannot extract transcript_id from $_";
            my $transcript_id = $1 or die "Error, no trans id from $_";
            
            my ($lend, $rend) = ($x[3], $x[4]);
            my $chr = $x[0];
            my $orient = $x[6];
            
            
            push (@{$genes{$chr}->{$gene_id}->{$transcript_id}}, { 
                gene => $gene_id,
                transcript => $transcript_id,
                chr => $chr,
                lend => $lend,
                rend => $rend,
                orient => $orient,
            }
                  );
            


        }
        close $fh;
    
        
        

    }

    
    my %target_to_aligns;

    open (my $fh, $gmap_gff3_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }

        chomp;
        my ($hit_acc, $filename, $type, $lend, $rend, $per_id, $orient, $dot, $info) = split(/\t/);

        my %info_hash;
        foreach my $keyval (split(/;/, $info)) {
            my ($key, $val) = split(/=/, $keyval);
            $info_hash{$key} = $val;
        }

        my $alignment_ID = $info_hash{ID};
        
        my ($target, $range_lend, $range_rend) = split(/\s+/, $info_hash{Target});

        push (@{$target_to_aligns{$target}->{$alignment_ID}}, { hit => $hit_acc,
                                                                lend => $lend,
                                                                rend => $rend,
                                                                orient => $orient,
                                                                per_id => $per_id,
                                                                
                                                                range_lend => $range_lend,
                                                                range_rend => $range_rend,
                                                            });
        
    }
    close $fh;

    
    ## header
    print join("\t", "#transcript", "num_alignments", "align_descr(s)", "[chim_annot_mapping]") . "\n";
    
    foreach my $target (keys %target_to_aligns) {
     
        ## create spans
        my @span_ids = keys %{$target_to_aligns{$target}};
        if (scalar @span_ids < 2) { next; } # only want candidate chimeras
        
        my @spans;
        foreach my $span_id (@span_ids) {
            my $exon_hits_aref = $target_to_aligns{$target}->{$span_id};
            my $span_struct = &convert_to_span($exon_hits_aref);
            push (@spans, $span_struct);
        }
        

        my $align_text = "$target";

        @spans = sort {$a->{range_lend}<=>$b->{range_lend}} @spans;
        
        my $num_aligns = scalar @spans;
        $align_text .= "\t$num_aligns";

        my $prev_align;
        
        my @chim_align_descrs;
        my @at_exon_junctions;
        
        foreach my $align (@spans) {
            my $hit = $align->{chr};
            
            
                    
            my $lend = $align->{lend};
            my $rend = $align->{rend};

            my $range_lend = $align->{range_lend};
            my $range_rend = $align->{range_rend};
            my $orient = $align->{orient};
            my $per_id = $align->{per_id};
            

            push (@chim_align_descrs, "[$hit:($range_lend-$range_rend)$lend-$rend ($orient) $per_id\%]");
            
            if ($prev_align) {
                if (%genes) {
                    my ($left_exon, $left_delta) = &examine_exon_junction("left", $prev_align);
                    my ($right_exon, $right_delta) = &examine_exon_junction("right", $align);
                    if ($left_exon && $right_exon) {
                        
                        my $left_exon_orient = $left_exon->{orient};
                        my $right_exon_orient = $right_exon->{orient};
                        
                        my $prev_align_orient = $prev_align->{orient};
                        my $curr_align_orient = $align->{orient};

                        ## both orient conflict, swap them (transcript is in antisense orientation)
                        if ($left_exon_orient ne $prev_align_orient && $right_exon_orient ne $curr_align_orient) {
                            
                            push (@at_exon_junctions, 
                                  $right_exon->{gene}, $right_delta,                                   
                                  $left_exon->{gene}, $left_delta, 
                                  join("--", 
                                       $right_exon->{gene},
                                       $left_exon->{gene}
                                       )
                                  );
                            
                        }
                        else {
                                                    
                            push (@at_exon_junctions, $left_exon->{gene}, $left_delta, $right_exon->{gene}, $right_delta, join("--", $left_exon->{gene}, $right_exon->{gene}));
                        }
                                                
                    }
                }
            }
            

            $prev_align = $align;

        }
        
        $align_text .= "\t" . join(";", @chim_align_descrs);
        
        unless (@at_exon_junctions) {
            @at_exon_junctions = ('.'); # placeholder
        }

        $align_text .= "\t" . join(";", @at_exon_junctions);
        
        
        print $align_text . "\n";
    }
    
    

    exit(0);
}

####
sub examine_exon_junction {
    my ($side_of_junction, $align_struct) = @_;

    #use Data::Dumper;
    #print Dumper($align_struct);
    #die;
    
    my $chr = $align_struct->{chr};
    my $lend = $align_struct->{lend};
    my $rend = $align_struct->{rend};
    my $orient = $align_struct->{orient};
    
    my $coord_to_examine; # lend or rend, depends on side_of_junction and orient
    if ($side_of_junction eq 'left') {
        if ($orient eq '+') {
            $coord_to_examine = "rend";
        }
        else {
            $coord_to_examine = "lend";
        }
    }
    elsif ($side_of_junction eq 'right') {
        if ($orient eq '+') {
            $coord_to_examine = 'lend';
        }
        else {
            $coord_to_examine = 'rend';
        }
    }

    my @hits;

    foreach my $gene_id (keys %{$genes{$chr}}) {
        
        foreach my $transcript_id (keys %{$genes{$chr}->{$gene_id}}) {
            
            my @exons = @{$genes{$chr}->{$gene_id}->{$transcript_id}};
            
            @exons = sort {$a->{lend}<=>$b->{lend}} @exons;
            
            ## exclude first and last exons, only looking at internal boundaries
            $exons[0]->{terminal} = 1;
            $exons[$#exons]->{terminal} = 1;
            
            my $num_exons = scalar(@exons);
            my $counter= 0;
            foreach my $exon (@exons) {
                $counter++;
                $exon->{exon_num} = "$counter/$num_exons";
            }
            
                        
            foreach my $exon (@exons) {
                if ($exon->{lend} < $rend && $exon->{rend} > $lend ) {

                    #push (@hits, "$gene_id:" . $exon->{lend} . "-" . $exon->{rend});
                    
                    my $delta = abs($align_struct->{$coord_to_examine} - $exon->{$coord_to_examine});
                    
                    #my $delta = &min( abs($align_struct->{lend} - $exon->{lend}),
                    #                  abs($align_struct->{rend} - $exon->{rend}),
                    #                  abs($align_struct->{lend} - $exon->{rend}),
                    #                 abs($align_struct->{rend} - $exon->{lend}),
                    #                );
                    
                    push (@hits, { delta => $delta,
                                   exon => $exon,

                                   # below for debugging
                                   pt_align => $align_struct->{$coord_to_examine},
                                   pt_exon => $exon->{$coord_to_examine},
                                   side_of_junction => $side_of_junction,
                                   orient => $orient,
                                   coord_to_examine => $coord_to_examine,
                                   
                               });
                    
                }
            }
            
        }
    }

    
    if (@hits) {

        #use Data::Dumper;  print Dumper(\@hits); print Dumper($align_struct);

        @hits = sort {$a->{delta}<=>$b->{delta}} @hits;
        my $top_hit = shift @hits;
        my $exon = $top_hit->{exon};
        my $delta = $top_hit->{delta};
        return($exon, $delta);
        #push (@hits, "$gene_id:" . $exon->{lend} . "-" . $exon->{rend});
    }
    else {
        return();
    }
    
    
        
}

####
sub min {
    my @vals = @_;

    @vals = sort {$a<=>$b} @vals;
    
    my $min_val = shift @vals;

    return($min_val);
}



####
sub convert_to_span {
    my ($exons_aref) = @_;

    my @chr_coords;
    my @span_coords;
    
    my $orient;
    my $chr;
    
    my $sum_per_id = 0;
    my $sum_len = 0;

    foreach my $exon (@$exons_aref) {
        push (@chr_coords, $exon->{lend});
        push (@chr_coords, $exon->{rend});
        push (@span_coords, $exon->{range_lend});
        push (@span_coords, $exon->{range_rend});

        my $len = abs($exon->{rend} - $exon->{lend} + 1);
        $sum_len += $len;
        $sum_per_id += $len * $exon->{per_id};
        

        my $exon_chr = $exon->{hit};
        if ($chr && $exon_chr ne $chr) {
            die "inconsistent chr assignments";
        }
        else {
            $chr = $exon_chr;
        }
        
        my $exon_orient = $exon->{orient};
        if ($orient && $exon_orient ne $orient) {
            die "inconsistent exon orient";
        }
        else {
            $orient = $exon_orient;
        }
    }
    
    @chr_coords = sort {$a<=>$b} @chr_coords;
    @span_coords = sort {$a<=>$b} @span_coords;
    
    my $chr_coords_lend = $chr_coords[0];
    my $chr_coords_rend = $chr_coords[$#chr_coords];
    
    my $span_coords_lend = $span_coords[0];
    my $span_coords_rend = $span_coords[$#span_coords];
    

    my $per_id = $sum_per_id / $sum_len;

    my $span_struct = { lend => $chr_coords_lend,
                        rend => $chr_coords_rend,
                        orient => $orient,
                        chr => $chr,

                        range_lend => $span_coords_lend,
                        range_rend => $span_coords_rend,

                        exons => $exons_aref,
                    
                        per_id => sprintf("%.2f", $per_id),
                        
                    };

    return($span_struct);
}


    
