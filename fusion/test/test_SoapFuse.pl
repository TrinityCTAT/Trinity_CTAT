#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;
use __GLOBALS__;

my $left_fq = "reads.left.simPE.fq.gz";
my $right_fq = "reads.right.simPE.fq.gz";

my $INSTALL_DIR = "$FindBin::Bin/../";

main: {

    
    my $pipeliner = new Pipeliner(-verbose => 1);
    
        
    
    ##############
    ## STAR-Fusion
    ##############
    
    $pipeliner->add_commands(
        new Command("$INSTALL_DIR/run_SoapFuse.pl reads.left.simPE.fq.gz reads.right.simPE.fq.gz soapfuse_outdir",
                    "soapfuse_outdir.ok")
        );
    
    
    ## Execute pipeline
    
    $pipeliner->run();

    
    
    exit(0);

}


                             



