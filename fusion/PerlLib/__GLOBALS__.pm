package __GLOBALS__;

use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT = qw($TRINITY_HOME $FUSION_ANNOTATOR_LIB);

our ($TRINITY_HOME, $FUSION_ANNOTATOR_LIB);



BEGIN {
    unless ($TRINITY_HOME = $ENV{TRINITY_HOME}) {
        confess "Error, need TRINITY_HOME env var set to Trinity base installation directory";
    }
    unless ($FUSION_ANNOTATOR_LIB = $ENV{FUSION_ANNOTATOR_LIB}) {
       confess "Error, need FUSION_ANNOTATOR_LIB env var set to fusion annotator lib resource dir";
    }
}


1;

#EOM

