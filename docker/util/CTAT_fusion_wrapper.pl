#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 left.fq.gz right.fq.gz ctat_genome_lib.tar.gz\n\n";

my $left_fq_gz = $ARGV[0] or die $usage;
my $right_fq_gz = $ARGV[1] or die $usage;
my $ctat_genome_lib_tar_gz = $ARGV[2] or die $usage;


my $STAR_FUSION_HOME = $ENV{STAR_FUSION_HOME} or die "Error, env var STAR_FUSION_HOME must be set";

main: {

    my $cmd = "tar xvf $ctat_genome_lib_tar_gz";
    &process_cmd($cmd);

    $cmd = "/usr/local/src/STAR-Fusion_v0.7.0_FULL/STAR-Fusion " .
	" --left_fq  $left_fq_gz" .
	" --right_fq $right_fq_gz" . 
        " --genome_lib_dir CTAT_lib" .
	" --output_dir star_fusion_outdir";

    &process_cmd($cmd);
    
    $cmd = "tar -zcvf star_fusion_outdir.tar.gz star_fusion_outdir";
    &process_cmd($cmd);

    print STDERR "Done.\n";
    
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
