#!/usr/bin/perl

### Imports ###################################################################
package ElementID;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

### Global Variables ##########################################################
my %options = (
    "help"   => 0,
    "man"    => 0,
    "bed"    => "",
    "inline" =>  ""
);

my $gencode_file     = "gencodev19.txt";
my $ctcf_file        = "ctcf_motifs.bed";
my $promoter_kb_dist = 3000;

### Subroutines ###############################################################
# parse_args
# Description:
#   parse command line input arguments
sub parse_args {
    GetOptions(
        'h|help'     => \$options{"help"},
        'man'        => \$options{"man"},
        'b|bed:s'    => \$options{"bed"},
        'i|inline:s' => \$options{"inline"}
    ) or pod2usage(2);
    if ($options{"help"}) {
        pod2usage(1);
    } elsif ($options{"man"}) {
        pod2usage(-exitval => 0, -verbose => 2);
    }
    # must have either bed or inline option
    pod2usage({ -message => q{Must have either '-b' or '-i' option},
                 -exitval => 1,
                 -verbose => 1
    }) unless ($options{"bed"} || $options{"inline"});
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
    my $readline = <$f_gencode>; # skip header line
    
    # test each line in gencode file
    while ($readline = <$f_gencode>) {
        chomp($readline);
        my @splitline = split(/\t/, $readline);
        my ($test_name, $test_chrom, $test_strand, $test_start, $test_end, $test_name2) = 
            @splitline;

        # if matching chromosomes do further work
        if ($chrom eq $test_chrom) {
            my $test_dist = distance($start, $end, $test_start, $test_end);
            # find the distance, and replace $dist if $test_dist is smaller
            if (!defined($dist) || $test_dist < $dist) {
                $closest_chrom = $test_chrom;
                $closest_start = $test_start;
                $closest_end   = $test_end;
                $dist          = $test_dist;
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

    # if intervals overlap
    if ($b_lower <= $a_upper && $b_lower >= $a_lower) {
        return(0);
    } elsif ($b_upper <= $a_upper && $b_upper >= $a_lower) {
        return(0);
    } else {
        # return smallest of the differences in bounds
        my $dist1 = abs($b_lower - $a_upper);
        my $dist2 = abs($a_lower - $b_upper);
        
        if ($dist1 < $dist2) {
            return($dist1);
        } else {
            return($dist2);
        }
    }
}

# is_promoter
# Description:
#   Check whether the given region is part of a promoter or not
# Inputs:
#   $chrom:    chromosome of region of interest
#   $start:    start site of region of interest
#   $end:      end site of region of interest
# Outputs:
#   True/False
sub is_promoter {
    my ($chrom, $start, $end) = @_;
    my ($tss_chrom, $tss_start, $tss_end, $tss_dist) = 
        closest_TSS($chrom, $start, $end);

    # if regions are detected within the window of being a promoter
    if ($tss_dist <= $promoter_kb_dist) {
        return(1);
    }
    return(0);
}

# is_anchor
# Description:
#   Check whether the given region is part of a chromatin loop anchor or not
# Inputs:
#   $chrom:    chromosome of region of interest
#   $start:    start site of region of interest
#   $end:      end site of region of interest
# Outputs:
#   True/False
sub is_anchor {
    my ($chrom, $start, $end) = @_;

    open(my $f_ctcf, "<", $ctcf_file) or die "Could not open $ctcf_file!\n";
    # test against each CTCF motif region
    while (my $readline = <$f_ctcf>) {
        chomp($readline);
        my @splitline = split(/\t/, $readline);
        my ($test_chrom, $test_start, $test_end, $test_strand, $test_fimo, $test_p, $test_q) = 
            @splitline;
        
        # if chromosomes are equal, and the distance between the regions is zero (i.e. the regions overlap)
        my $dist = distance($start, $end, $test_start, $test_end);

        if ($chrom eq $test_chrom && 0 == $dist) {
            close($f_ctcf);
            return(1);
        }
    }
    close($f_ctcf);
    return(0);
}

### Main ######################################################################
parse_args();

if ($options{"inline"}) {
    my @region_list = split(/:/, $options{"inline"});
    foreach my $region (@region_list) {
        my ($chrom, $start, $end) = split(/,/, $region);
        my $is_prom = is_promoter($chrom, $start, $end);
        my $is_anch = is_anchor($chrom, $start, $end);  
        print(join("\t", $chrom, $start, $end, $is_prom, $is_anch), "\n");
    }
} else {
    open(my $f_in, "<", $options{"bed"}) or die "Could not open $options{'bed'}!\n";
    while (my $readline = <$f_in>) {
        chomp($readline);
        # skip header line if it exists
        if (grep(/^#/, $readline) || grep(/^chrom/, $readline)) {
            next;
        }

        my @splitline = split(/\t/, $readline);
        my ($chrom, $start, $end) = @splitline;
        my $is_prom = is_promoter($chrom, $start, $end);
        my $is_anch = is_anchor($chrom, $start, $end);  
        print(join("\t", $chrom, $start, $end, $is_prom, $is_anch), "\n");
    }
    close($f_in);
}

__END__

=head1 NAME

elementid: Description

=head1 SYNOPSIS

perl elementid.pl -b <BED file>

perl elementid.pl -i <chrom>,<start>,<end>[:<chrom>,<start>,<end>...]

Options:
    -h, --help          Brief help message
    --man               Man page with full documentation
    -b, --bed           BED file containing regions to check
    -i, --inline        Check region from commandline. Chromosome, start, and end should be comma-separated.
                        Multiple regions should be colon-separated.