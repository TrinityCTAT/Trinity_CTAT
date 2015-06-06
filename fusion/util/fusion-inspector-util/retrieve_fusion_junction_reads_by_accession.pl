#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Data::Dumper;
use Carp;

my $usage = "usage: $0 read_names.accs  fileA.bam,fileB.bam,...\n\n";

my $read_name_accs_list = $ARGV[0] or die $usage;
my $bam_file_listing = $ARGV[1] or die $usage;

main: {
    
    my %reads_want;
    {
        my @read_infos = `cat $read_name_accs_list | egrep -v ^\#`;
        chomp @read_infos;
        foreach my $read_info (@read_infos) {
            my ($geneA, $geneB, $reads_list) = split(/\t/, $read_info);
            my $fusion_contig = "$geneA--$geneB";
            foreach my $read_name (split(/,/, $reads_list)) {
                $reads_want{"$fusion_contig|$read_name"} = 1;
            }
        }
                
    }

    foreach my $bam_file (split(/,/, $bam_file_listing) ) {
        
        my $sam_reader = new SAM_reader($bam_file);
        while (my $sam_entry = $sam_reader->get_next()) {
            my $scaffold = $sam_entry->get_scaffold_name();
            my $read_name = $sam_entry->reconstruct_full_read_name();
            $read_name = "$scaffold|$read_name";
            if ($reads_want{$read_name}) {
                print $sam_entry->get_original_line() . "\n";
                delete($reads_want{$read_name});
            }
        }

    }

    if (%reads_want) {
        confess "Error, missing junction reads: " . Dumper(\%reads_want);
        
    }
    
    exit(0);
}


