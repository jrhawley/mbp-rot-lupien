#!/usr/bin/perl

### Imports ###################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

### Global Variables ##########################################################
my %options = ();

### Subroutines ###############################################################
# parse_args
# Description:
#   parse command line input arguments
sub parse_args {
    GetOptions(
        \%options,
        'h|help',
        'man'
    ) or pod2usage(2);
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage(-exitval => 0, -verbose => 2);
    }
}

### Main ######################################################################
parse_args();

__END__

=head1 NAME

elementid: Description

=head1 SYNOPSIS

perl elementid.pl -b <BED file>

Options:
    -h, --help          Brief help message
    --man               Man page with full documentation
    -b                  BED file containing regions to check