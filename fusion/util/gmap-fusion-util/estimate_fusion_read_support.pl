#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use File::Basename;
use Process_cmd;
use Pipeliner;
use SAM_reader;

my $usage = "\n\n\tusage: $0 trans.fasta gmap.map.gff3.chims_described left.fq right.fq\n\n";

my $trans_fasta = $ARGV[0] or die $usage;
my $chims_described = $ARGV[1] or die $usage;
my $left_fq_file = $ARGV[2] or die $usage;
my $right_fq_file = $ARGV[3] or die $usage;

my $ANCHOR = 12;


$trans_fasta = &ensure_full_path($trans_fasta);
$chims_described = &ensure_full_path($chims_described);
$left_fq_file = &ensure_full_path($left_fq_file);
$right_fq_file = &ensure_full_path($right_fq_file);


foreach my $file ($trans_fasta, $chims_described, $left_fq_file, $right_fq_file) {
    unless (-s $file) {
        confess "Error, cannot locate file $file";
    }
}


main: {

    my %chims = &parse_chims($chims_described);
    
    my $chim_candidates_fasta_filename = basename($chims_described) . ".fasta";
    
    &extract_chim_candidate_seqs($trans_fasta, $chim_candidates_fasta_filename, \%chims);
    
    my $bowtie2_bam = &align_reads_using_bowtie2($chim_candidates_fasta_filename, $left_fq_file, $right_fq_file);
        
    my %fusion_support = &capture_fusion_support($bowtie2_bam, \%chims);

    
    ## generate output, include junction and spanning frag support info:
    
    foreach my $target_trans_id (keys %chims) {
        my $chim_info = $chims{$target_trans_id};
        my $line = $chim_info->{line};

        my @junction_reads;
        my @spanning_frags;
        
        if (exists $fusion_support{$target_trans_id}) {
            my $info_href = $fusion_support{$target_trans_id};
            if (exists $info_href->{spanning}) {
                @spanning_frags = keys %{$info_href->{spanning}};
            }
            if (exists $info_href->{junction}) {
                @junction_reads = keys %{$info_href->{junction}};
            }

        }
        
        my $spanning_frag_list = join(",", @spanning_frags) || ".";
        my $junction_frag_list = join(",", @junction_reads) || ".";

        my $J = scalar(@junction_reads);
        my $S = scalar(@spanning_frags);

        print join("\t", $line, $J, $S, $junction_frag_list, $spanning_frag_list) . "\n";
    }
    
    exit(0);

}


####
sub capture_fusion_support {
    my ($bam_file, $chims_href) = @_;


    my %fusion_support;

    # sort according to read name and contig
    # these should all be perfect pairs

    my $pipeliner = new Pipeliner(-verbose => 2);

    my $sorted_alignments_file = "$bam_file.sorted_by_read_n_contig.sam.gz";
    my $cmd = "samtools view $bam_file | sort -k1,1 -k3,3 | gzip -c - > $sorted_alignments_file";

    $pipeliner->add_commands(new Command($cmd, "$sorted_alignments_file.ok"));
    
    $pipeliner->run();

    
    my $sam_reader = new SAM_reader($sorted_alignments_file);

    while ($sam_reader->has_next()) {
        my $sam_entryA = $sam_reader->get_next();
        my $sam_entryB = $sam_reader->get_next();
        
        my $target_trans_id = $sam_entryA->get_scaffold_name();
        unless ($sam_entryB->get_scaffold_name() eq $target_trans_id) {
            confess "Error, trans IDs dont match for supposed pairs.";
        }

        my $frag_name = $sam_entryA->get_core_read_name();
        if ($frag_name ne $sam_entryB->get_core_read_name()) {
            confess "Error, core frag names dont match up: $frag_name vs. " . $sam_entryB->get_core_read_name();
        }
        

        unless ($sam_entryA->is_first_in_pair() xor $sam_entryB->is_first_in_pair()) {
            print STDERR "error, have unpaired pair... skipping.\n";
            next;
        }
        
        
        my $brkpt_range = $chims_href->{$target_trans_id}->{brkpt_range};
        my ($break_left, $break_right) = split(/-/, $brkpt_range);


        my ($trans_coords_A_aref, @trash1) = $sam_entryA->get_alignment_coords();
        my ($trans_coords_B_aref, @trash2) = $sam_entryB->get_alignment_coords();

        # sort them according by coordinate
        ($trans_coords_A_aref, $trans_coords_B_aref) = sort {$a->[0]->[0] <=> $b->[0]->[0]} ($trans_coords_A_aref, $trans_coords_B_aref);

        my ($A_lend, $A_rend) = ($trans_coords_A_aref->[0]->[0], $trans_coords_A_aref->[$#$trans_coords_A_aref]->[1]);
        my ($B_lend, $B_rend) = ($trans_coords_B_aref->[0]->[0], $trans_coords_B_aref->[$#$trans_coords_B_aref]->[1]);
        
        if ($A_lend < $break_left && $B_rend > $break_right) {
            ## fragment overlaps breakpoint.

            ## determine if it's a spanning pair or a fusion junction.
            if ($A_rend < $break_left && $B_lend > $break_right) {
                ## a spanning pair:

                #print STDERR "Found spanning fragment for $target_trans_id\n";
                $fusion_support{$target_trans_id}->{spanning}->{$frag_name}++;
            }
            else {
                ## see if any alignment overlaps the point of the junction
                foreach my $align_seg (@$trans_coords_A_aref, @$trans_coords_B_aref) {
                    my ($lend, $rend) = @$align_seg;
                    
                    ## ensure the alignment meets the anchor requirement.

                    #                     brktp
                    #    ------------------|------------------------ 
                    #           <-- anchor on each side --->
                    
                    
                    if ($lend <= ($break_left - $ANCHOR) && $rend >= ($break_right + $ANCHOR)) {
                        # overlaps junction breakpoint
                        #print STDERR "Found a JUNCTION read for $target_trans_id\n";
                        $fusion_support{$target_trans_id}->{junction}->{$frag_name}++;
                        last;
                    }
                }
            }
        }
    }
    
    return(%fusion_support);

}



####
sub align_reads_using_bowtie2 {
    my ($trans_fasta, $left_fq_file, $right_fq_file) = @_;

    my $cmd = "ln -sf $trans_fasta bowtie2_target.fa";
    &process_cmd($cmd);

    my $pipeliner = new Pipeliner(-verbose=>1);
    $cmd = "bowtie2-build bowtie2_target.fa bowtie2_target > /dev/null";
    $pipeliner->add_commands(new Command($cmd, "bowtie2_target.build.ok"));

    $cmd = "bash -c \"set pipefail -o && bowtie2 -k10 -p 4 --no-mixed --no-discordant --very-fast --local -x bowtie2_target -1 $left_fq_file -2 $right_fq_file "
        . " | samtools view -F 4 -Sb - | samtools sort -@ 4 - $trans_fasta.bowtie2\"";
    $pipeliner->add_commands(new Command($cmd, "bowtie2_align.ok"));
    
    $pipeliner->run();
    
    return("$trans_fasta.bowtie2.bam");
    
}



####
sub extract_chim_candidate_seqs {
    my ($trans_fasta_filename, $output_filename, $chims_href) = @_;

    my $checkpoint = "$output_filename.ok";
    if (-e $checkpoint && -s $trans_fasta_filename) {
        return;
    }
    else {
        my $fasta_reader = new Fasta_reader($trans_fasta_filename);
        open (my $ofh, ">$output_filename") or die "Error, cannot write to file: $trans_fasta_filename";
        
        while (my $seq_obj = $fasta_reader->next()) {
            
            my $accession = $seq_obj->get_accession();
            
            if (exists $chims_href->{$accession}) {
                
                my $sequence = $seq_obj->get_sequence();
                
                print $ofh ">$accession\n$sequence\n";
            }
        }
        
        close $ofh;
        
        &process_cmd("touch $checkpoint");
        
        return;
    }
}
        



####
sub parse_chims {
    my ($chims_described_file) = @_;

    my %chims;
    
    open (my $fh, $chims_described_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; } # header or comment
        chomp;
        my $line = $_;

        my @x = split(/\t/);
        
        my $trans_acc = $x[0];
        my $fusion_info = $x[3];

        my ($geneA, $deltaA, $trans_brkptA, 
            $chrA_n_coordA,
            $geneB, $deltaB, $trans_brkptB, 
            $chrB_n_coordB,
            $fusion_name) = split(/;/, $fusion_info);

        my $brkpt_range = join("-", sort ($trans_brkptA, $trans_brkptB));
        
        $chims{$trans_acc} = { line => $line,
                               brkpt_range => $brkpt_range,
        };
    }
    close $fh;

    return(%chims);
}


                               
            
          
