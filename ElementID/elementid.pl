#!/usr/bin/perl

### Imports ###################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

### Global Variables ##########################################################
my %options = ();

my $gencode_file = "gencodev19.txt";

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

# closest_TSS
# Description:
#   Given a chromosome and position, find the closest transcription start site for a gene
# Inputs:
#   chrom:              chromosome
#   start:              region start site
#   end:                region end site
#   strand:             strand of testing region
# Outputs:
#   closest_chrom:      chromosome of closest TSS
#   closest_start:      start site of closest TSS
#   closest_end:        end site of closest TSS
#   closest_strand:     strand of closest TSS
#   distance:           distance between region of interest and closest TSS
sub closest_TSS {
    my ($chrom, $start, $end) = @_;
    my ($closest_chrom, $closest_start, $closest_end, $dist);
    open(my $f_gencode, "<", $gencode_file) or die "Could not open $gencode_file!\n";
    my $readline = <$f_gencode>;
    while ($readline = <$f_gencode>) {
        chomp($readline);
        my @splitline = split(/\t/, $readline);
        my $test_chrom = $splitline[1];
        my $test_start = $splitline[3];
        my $test_end   = $splitline[4];

        if ($chrom eq $test_chrom) {
            my $test_dist = distance($start, $end, $test_start, $test_end);
            if ($test_dist > $dist) {
                $dist = $test_dist;
            }
        }
    }
    close($f_gencode);

    return($closest_chrom, $closest_start, $closest_end, $dist);
}

# distance
# Description:
#   Find shortest distance between two intervals. Returns 0 if intervals overlap
# Inputs:
#   a_lower:    lower bound of interval a
#   a_upper:    upper bound of interval a
#   b_lower:    lower bound of interval b
#   b_upper:    upper bound of interval b
# Outputs:
#   dist:       shortest distance between intervals a and b
sub distance {
    my ($a_lower, $a_upper, $b_lower, $b_upper) = @_;

    if ($b_lower <= $a_upper && $b_lower >= $a_lower) {
        return(0);
    } elsif ($b_upper <= $a_upper && $b_upper >= $a_lower) {
        return(0);
    } else {
        my $dist1 = abs($b_lower - $a_upper);
        my $dist2 = abs($a_lower - $b_upper);
        
        if ($dist1 > $dist2) {
            return($dist1);
        } else {
            return($dist2);
        }
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