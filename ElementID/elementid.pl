#!/usr/bin/perl

### Imports ###################################################################
package ElementID;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

### Global Variables ##########################################################
my %options = (
    "help" => 0,
    "man"  => 0,
    "bed"  => ""
);

my $gencode_file = "gencodev19.txt";

### Subroutines ###############################################################
# parse_args
# Description:
#   parse command line input arguments
sub parse_args {
    GetOptions(
        'h|help' => \$options{"help"},
        'man'    => \$options{"man"},
        'b|bed=s'  => \$options{"bed"}
    ) or pod2usage(2);
    if ($options{"help"}) {
        pod2usage(1);
    }
    elsif ($options{"man"}) {
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
# Outputs:
#   closest_chrom:      chromosome of closest TSS
#   closest_start:      start site of closest TSS
#   closest_end:        end site of closest TSS
#   distance:           distance between region of interest and closest TSS
sub closest_TSS {
    my ($chrom, $start, $end) = @_;
    my ($closest_chrom, $closest_start, $closest_end, $dist);
    open(my $f_gencode, "<", $gencode_file) or die "Could not open $gencode_file!\n";
    my $readline = <$f_gencode>;
    while ($readline = <$f_gencode>) {
        chomp($readline);
        my @splitline = split(/\t/, $readline);
        my ($test_name, $test_chrom, $test_strand, $test_start, $test_end, $test_name2) = 
            @splitline;

        # print($test_chrom, "\n");
        if ($chrom eq $test_chrom) {
            # print("[$start, $end]\t[$test_start, $test_end]\n");
            my $test_dist = distance($start, $end, $test_start, $test_end);
            # print($test_dist, "\n");
            if (!defined($dist) || $test_dist < $dist) {
                $closest_chrom = $test_chrom;
                $closest_start = $test_start;
                $closest_end   = $test_end;
                $dist          = $test_dist;
            }
        }
        # sleep(1);
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
        
        if ($dist1 < $dist2) {
            return($dist1);
        } else {
            return($dist2);
        }
    }
}

### Main ######################################################################
unless (caller) {
    parse_args();
    open(my $f_in, "<", $options{"bed"}) or die "Could not open $options{'bed'}!\n";
    while (my $readline = <$f_in>) {
        chomp($readline);
        if (grep(/^#/, $readline) || grep(/^chrom/, $readline)) {
            next;
        }

        my @splitline = split(/\t/, $readline);
        my ($chrom, $start, $end) = @splitline;
        print(join("\t", closest_TSS($chrom, $start, $end)), "\n");
    }
    close($f_in);
}

__END__

=head1 NAME

elementid: Description

=head1 SYNOPSIS

perl elementid.pl -b <BED file>

Options:
    -h, --help          Brief help message
    --man               Man page with full documentation
    -b                  BED file containing regions to check