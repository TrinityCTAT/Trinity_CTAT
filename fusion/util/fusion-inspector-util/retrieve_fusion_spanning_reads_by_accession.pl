#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;

my $usage = "usage: $0 read_names.accs  fileA.bam,fileB.bam,...\n\n";

my $read_name_accs_list = $ARGV[0] or die $usage;
my $bam_file_listing = $ARGV[1] or die $usage;

main: {
    
    my %cores_want;
    {
        my @read_infos = `cat $read_name_accs_list | egrep -v ^\#`;
        chomp @read_infos;
        foreach my $read_info (@read_infos) {
            my ($geneA, $geneB, $reads_list) = split(/\t/, $read_info);
            my $fusion_contig = "$geneA--$geneB";
            foreach my $read_name (split(/,/, $reads_list)) {
                $cores_want{"$fusion_contig|$read_name"} = 1;
            }
        }
                
    }

    
    my %reads_seen;

    foreach my $bam_file (split(/,/, $bam_file_listing) ) {
        
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            my $scaffold = $sam_entry->get_scaffold_name();
            
            my $core_read_name = $sam_entry->get_core_read_name();
            
            $core_read_name = "$scaffold|$core_read_name";
            if ($cores_want{$core_read_name}) {
                
                
                my $full_read_name = $sam_entry->reconstruct_full_read_name();
                $full_read_name = "$scaffold|$full_read_name";
                $full_read_name =~ m|/([12])$| or die "Error cannot parse read end from $full_read_name";
                my $end = $1;
                my $opposite_end = ($end == 1) ? 2 : 1;
                my $opposite_read_name = $core_read_name . "/" . $opposite_end; 
                if (! $reads_seen{$full_read_name}) {
                    $reads_seen{$full_read_name}++;
                    print $sam_entry->get_original_line() . "\n";
                
                    if ($reads_seen{$opposite_read_name}) {
                        delete $cores_want{$core_read_name};
                    }
                }
                
            }
        }

    }

    if (%cores_want) {
        die "Error, missing spanning reads: " . Dumper(\%cores_want);
        
    }
    
    exit(0);
}


