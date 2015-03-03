#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;
use __GLOBALS__;

my $left_fq = "reads.left.simPE.fq.gz";
my $right_fq = "reads.right.simPE.fq.gz";

main: {

    
    my $pipeliner = new Pipeliner(-verbose => 1);
    
    ####################################
    ## run Trinity, assemble transcripts
    ####################################

    my $cmd = "$TRINITY_HOME/Trinity --left $left_fq --right $right_fq --seqType fq "
        . " --max_memory 2G --CPU 2 --output trinity_out_dir --verbose --full_cleanup ";

    $pipeliner->add_commands(new Command($cmd, "trinity_out_dir.ok"));

    



    $pipeliner->run();

    
    exit(0);
}


                             



