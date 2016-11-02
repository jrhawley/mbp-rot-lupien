#!/usr/bin/perl

### Imports ###################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

### Global Variables ##########################################################
my $man              = 0;
my $help             = 0;
my $c3d_file         = "";
my $bed_file         = "";
my $output_directory = "";
my $window           = 0;

### Subroutines ###############################################################
# parse_args
# Description:
#   parse command line input arguments
sub parse_args {
    GetOptions(
        'b|bed=s'    => \$bed_file,
        'c|c3d=s'    => \$c3d_file,
        'o|output=s' => \$output_directory,
        'w|window=i' => \$window,
        'help|?'     => \$help,
        'man'        => \$man
    ) or pod2usage(2);
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage(-exitval => 0, -verbose => 2);
    }
    pod2usage({ -message => q{Mandatory argument '-b' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $bed_file;
    pod2usage({ -message => q{Mandatory argument '-c' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $c3d_file;
    pod2usage({ -message => q{Mandatory argument '-w' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $window;
}

### Main ######################################################################
parse_args();

# store gene list info in hash
print("Reading genes BED file\n");
my %genes = ();
open(my $f_bed, "<", $bed_file) or die "Could not open $bed_file\n";
while (<$f_bed>) {
    chomp();
    my @line = split(/\t/);
    $genes{$line[0]} = {
        "chrom" => $line[1],                        # chromosome
        "start" => $line[2],                        # start position
        "end"   => $line[3],                        # end position
        "path"  => join(".", $c3d_file, $line[0]),  # output file containing filtered C3D output 
        "io"    => undef                            # output file I/O handler
    };
}
close($f_bed);

# sort through C3D file and extract lines of interest
# I'm assuming that the C3D file will be much larger than the genes file, so I'm structuring the loops this way to be efficient
print("Reading C3D file\n");
foreach my $gene (keys(%genes)) {
    open($genes{$gene}{"io"}, ">", $genes{$gene}{"path"});
}
open(my $f_c3d, "<", $c3d_file) or die "Could not open $c3d_file\n";
while (<$f_c3d>) {
    my $readline = $_;
    chomp();
    my @line = split(/\t/);
    my @chrom_region = split(/[:-]/, $line[0]);
    foreach my $gene (keys(%genes)) {
        # check if region in C3D file is within the local window of a gene of interest
        if ($chrom_region[0] eq $genes{$gene}{"chrom"} &&
            # start between gene boundaries
            ((($chrom_region[1] > $genes{$gene}{"start"} - $window) && ($chrom_region[1] < $genes{$gene}{"end"} + $window)) ||
            # end between gene boundaries
            (($chrom_region[2] > $genes{$gene}{"start"} - $window) && ($chrom_region[2] < $genes{$gene}{"end"} + $window)) ||
            # region encompasses gene boundaries
            (($chrom_region[1] < $genes{$gene}{"start"} - $window) && ($chrom_region[2] > $genes{$gene}{"end"} + $window))
        )) {
            print ${$genes{$gene}{"io"}} $readline;
        }
    }
}
close($f_c3d);
foreach my $gene (keys(%genes)) {
    close($genes{$gene}{"io"});
}

__END__

=head1 NAME

split-c3d: grab relevant C3D output lines based on gene(s) of interest, and window

=head1 SYNOPSIS

perl split-c3d.pl -b <gene(s) BED> -c <C3D output file> -w <window> [options]

Options:
    --help              Brief help message
    --man               Man page with full documentation
    -o, --output        Directory path for output files