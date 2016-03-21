#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 left.fq.gz right.fq.gz ctat_genome_lib.tar.gz\n\n";

my $left_fq_gz = $ARGV[0] or die $usage;
my $right_fq_gz = $ARGV[1] or die $usage;
my $ctat_genome_lib_tar_gz = $ARGV[2] or die $usage;


my $STAR_FUSION_HOME = $ENV{STAR_FUSION_HOME} or die "Error, env var STAR_FUSION_HOME must be set";
my $FUSION_INSPECTOR_HOME = $ENV{FUSION_INSPECTOR_HOME} or die "Error, env var FUSION_INSPECTOR_HOME must be set";


main: {

    my $cmd = "tar xvf $ctat_genome_lib_tar_gz";
    &process_cmd($cmd) unless (-d "CTAT_lib");
    
    ## Run STAR-Fusion
    
    $cmd = "$STAR_FUSION_HOME/STAR-Fusion " .
	" --left_fq  $left_fq_gz" .
	" --right_fq $right_fq_gz" . 
        " --genome_lib_dir CTAT_lib" .
	" --output_dir star_fusion_outdir";

    &process_cmd($cmd) unless (-d "star_fusion_outdir");

    ## Run FusionInspector
    $cmd = "$FUSION_INSPECTOR_HOME/FusionInspector --fusions star_fusion_outdir/star-fusion.fusion_candidates.final.abridged.FFPM " .
        " --genome_lib CTAT_lib " .
        " --left_fq $left_fq_gz " .
        " --right $right_fq_gz " .
        " --out_dir Fusion_Inspector-STAR " .
        " --out_prefix finspector " .
        " --align_utils STAR --prep_for_IGV --no_cleanup ";
    
    &process_cmd($cmd) unless (-d "Fusion_Inspector-STAR");
    
    
    ## Save outputs:
    my $ctat_outdir = "ctat_out";
    unless (-d $ctat_outdir) {
        mkdir "$ctat_outdir" or die "Error, cannot mkdir $ctat_outdir";
    }

    $cmd = "cp star_fusion_outdir/star-fusion.fusion_candidates.final.abridged.FFPM $ctat_outdir";
    &process_cmd($cmd);

    $cmd = "cp Fusion_Inspector-STAR/finspector.fusion_predictions.final.abridged.FFPM $ctat_outdir";
    &process_cmd($cmd);

    $cmd = "tar -zcvf $ctat_outdir.tar.gz $ctat_outdir";
    &process_cmd($cmd);
    
    exit(0);
    
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);
    if ($ret) {
	die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}
