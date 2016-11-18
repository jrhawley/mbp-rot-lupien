#!/usr/bin/perl

### Imports ###################################################################
package MuSE::MuSE;
use strict;
use warnings;
use Statistics::R;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path remove_tree);
use File::Basename;
use Parallel::ForkManager;
use Cwd;


### Global Variables ##########################################################
# Runtime options
my %options = (
    "Help"             => 0,
    "Manual"           => 0,
    "Mutations File"   => "",
    "C3D File"         => "",
    "DHS File"         => "",
    "Window"           => 0,    # Window size '-w' must be equal or smaller than the C3D window size used
    "Threshold"        => 0,
    "Output Directory" => getcwd(),
    "Parallel"         => 1
);

my $max_lines_split  = 10000;       # number of lines for a single C3D split block
my $max_processes    = 12;
my $prefix_muse      = "tmpmuse";   # prefix for temporary MuSE files
my $suffix_muse      = "MUSE";
my $suffix_mut  = "MUT";

my @tmp_filelist     = ();
my $header_line = join("\t",
        "GENE",
        "M_TYPE_1",
        "M_TYPE_1_BMR",
        "M_TYPE_1_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_1_N_MUT",
        "M_TYPE_1_P",
        "M_TYPE_1_P_LOCAL",
        "M_TYPE_2",
        "M_TYPE_2_BMR",
        "M_TYPE_2_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_2_N_MUT",
        "M_TYPE_2_P",
        "M_TYPE_2_P_LOCAL",
        "M_TYPE_3",
        "M_TYPE_3_BMR",
        "M_TYPE_3_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_3_N_MUT",
        "M_TYPE_3_P",
        "M_TYPE_3_P_LOCAL",
        "M_TYPE_4",
        "M_TYPE_4_BMR",
        "M_TYPE_4_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_4_N_MUT",
        "M_TYPE_4_P",
        "M_TYPE_4_P_LOCAL",
        "M_TYPE_5",
        "M_TYPE_5_BMR",
        "M_TYPE_5_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_5_N_MUT",
        "M_TYPE_5_P",
        "M_TYPE_5_P_LOCAL",
        "M_TYPE_6",
        "M_TYPE_6_BMR",
        "M_TYPE_6_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_6_N_MUT",
        "M_TYPE_6_P",
        "M_TYPE_6_P_LOCAL",
        "M_TYPE_7",
        "M_TYPE_7_BMR",
        "M_TYPE_7_BMR_LOCAL",
        "SIME_LENGTH",
        "M_TYPE_7_N_MUT",
        "M_TYPE_7_P",
        "M_TYPE_7_P_LOCAL",
        "FISHER_P",
        "FISHER_P_LOCAL",
        "N_DHS"
    );

### Subroutines ###############################################################
# binomial
# Description:
#   perform binomial test using R
# Inputs:
#   x: reference allele count
#   n: total minus errors
#   p: observe or expected frequency in genomic DNA
# Outputs:
#   p_value: p-value resulting from binomial test
sub binomial {
    my $x = shift;
    my $n = shift;
    my $p = shift;

    # start R instance
    my $R = Statistics::R->new();

    # set variables
    $R->set('x', $x);
    $R->set('n', $n);
    $R->set('p', $p);

    # perform test
    $R->run(q`b <- binom.test(x, n, p, alternative="greater")`);
    my $p_value = $R->get('b$p.value');

    $R->stop();

    return($p_value);
}

# fisher_combine_pval
# Description:
#   pvals: list of p-values
# Inputs:
#   pvals: array of p-values to be combined
# Outputs:
#   p_value: combined p-value
sub fisher_combine_pval {
    #_COMMENT_: why not initially cast as array?
    my $pvals = shift;
    my $R = Statistics::R->new();

    # combine values
    $R->set('pvals', [@{$pvals}]);
    $R->run(q`chisq <- -2*sum(log(pvals))`);
    $R->run(q`pval.combined <- 1 - pchisq(chisq, df = 2*length(pvals))`);

    my $p_value = $R->get('pval.combined');

    return($p_value);
}

# log10
# Description:
#   calculate log base 10
# Inputs:
#   n: input value
# Outputs:
#   log_10(n)
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

# assign_type
# Description:
#   Assign a generic type of mutation given the reference and mutant alleles
# Inputs:
#   a: reference allele
#   b: mutant allele
# Outputs:
#   type: type of mutation
sub assign_type {
    my $a = shift;
    my $b = shift;
    my $type;

    if (($a eq "A" && $b eq "C") || ($a eq "T" && $b eq "G")) {
        $type = "ACTG";
    } elsif (($a eq "A" && $b eq "G") || ($a eq "T" && $b eq "C")) {
        $type="AGTC";
    } elsif (($a eq "A" && $b eq "T") || ($a eq "T" && $b eq "A")) {
        $type="ATTA";
    } elsif (($a eq "C" && $b eq "A") || ($a eq "G" && $b eq "T")) {
        $type="CAGT";
    } elsif (($a eq "C" && $b eq "G") || ($a eq "G" && $b eq "C")) {
        $type="CGGC";
    } elsif (($a eq "C" && $b eq "T") || ($a eq "G" && $b eq "A")) {
        $type="CTGA";
    } else {
        $type="OTHER";
    }
    #_COMMENT_: add "NONE" mutation? just for robustness
    return($type);
}

# parse_args
# Description:
#   parse command line input arguments
sub parse_args {
    GetOptions(
        \%options,
        'help|h'           => \$options{"Help"},
        'man'              => \$options{"Manual"},
        'm|mutations=s'    => \$options{"Mutations File"},
        'c|c3d=s'          => \$options{"C3D File"},
        'd|dhs=s'          => \$options{"DHS File"},
        'w|window=i'       => \$options{"Window"},
        't|threshold=f'    => \$options{"Threshold"},
        'o|output:s'       => \$options{"Output Directory"},
        'parallel!'        => \$options{"Parallel"}
    ) or pod2usage(2);

    if ($options{"Help"}) {
        pod2usage(1);
    }
    elsif ($options{"Manual"}) {
        pod2usage(-exitval => 0, -verbose => 2);
    }
    pod2usage({ -message => q{Mandatory argument '-m' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $options{"Mutations File"};
    pod2usage({ -message => q{Mandatory argument '-c' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $options{"C3D File"};
    pod2usage({ -message => q{Mandatory argument '-b' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $options{"DHS File"};
    pod2usage({ -message => q{Mandatory argument '-w' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $options{"Window"};
    pod2usage({ -message => q{Mandatory argument '-t' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $options{"Threshold"};

    if (!(-e $options{"Output Directory"} && -d $options{"Output Directory"})) {
        make_path($options{"Output Directory"});
    }
}

# parse_mutations
# Description:
#   Read mutations file and store in hash structure
# Inputs:
#   mutation_file:  mutation file to be read
# Outputs:
#   mutations:      pointer to array for each line in mutations file
#   starts:         pointer to array containing line numbers of starting positions
#   chroms:         pointer to array containing chromosomes in use
sub parse_mutations {
    my $mutation_file = shift;

    my @mutations;
    my @starts;
    my @chroms;
    my %landmarks = ();
    my $n         = -1; # array position
    my %mutP      = ();

    open(my $inputMUT, "<", $mutation_file) or die "Could not open $mutation_file!\n";
    while (<$inputMUT>) {
        chomp();
        $n += 1;
        push(@mutations, $_);

        # extract information from mutations file
        my @line = split(/\t/);
        my $chr  = $line[0];
        my $bp   = $line[2];
        my $id   = $line[3];
        my $A1   = $line[4];
        my $A2   = $line[5];

        # save locations of points of interest, and the chromosomes their on
        if (!exists $landmarks{$chr}) {
            $landmarks{$chr} = $n;
            push(@starts, $n);
            push(@chroms, $chr);
        }  

        # add to mutation count for that sample
        if (!exists $mutP{$id}) {
            $mutP{$id} = 1;
        } else {
            $mutP{$id} += 1;
        }
    }
    close($inputMUT);

    return(\@mutations, \@starts, \@chroms);
}

# parse_reference
# Description:
#   Read reference BED file and store information
# Inputs:
#   reference_file:     reference BED file path
#   chroms_ref:         pointer to @chroms
#   starts_ref:         pointer to @starts 
#   mutations_ref:      pointer to @mutations
# Outputs:
#   wgdist:             distance of genome analyzed in reference BED (for calculating global mutation rates)
#   wgbmr:              pointer to hash containing whole genome background mutation rates by mutation type
sub parse_reference {
    my $reference_file = shift;
    my $chroms_ref     = shift;
    my @chroms         = @$chroms_ref;
    my $starts_ref     = shift;
    my @starts         = @$starts_ref;
    my $mutations_ref  = shift;
    my @mutations      = @$mutations_ref;

    my $wgdist  = 0;
    my %wgbmr = (
        "ACTG"  => 0,
        "AGTC"  => 0,
        "ATTA"  => 0,
        "CAGT"  => 0,
        "CGGC"  => 0,
        "CTGA"  => 0,
        "OTHER" => 0,
    );

    open(my $inputBED, "<", $reference_file) or die "Could not open $reference_file!\n";
    #_EXPLANATION_
    while (<$inputBED>) {
        chomp();
        my @line = split(/\t/);
        my $dist = $line[2] - $line[1];
        $wgdist += $dist;

        my $start = 0;
        my $fs    = 0;
        my $fn    = 0;

        for (my $i = 0; $i <= $#chroms; $i++) {
            if ($line[0] eq $chroms[$i]) {
                if ($i != $#chroms) {
                    $fs = $starts[$i];
                    $fn = $starts[$i + 1];
                } else {
                    $fs = $starts[$i];
                    $fn = $#mutations;
                }
            }
        }

        my $d = 0;
        my $t = 0;
        my $s = 10;
        while ( $d < $s ) {
            #_BAD_NAME_
            my @ps = split(/\t/, $mutations[$fs]);
            #_BAD_NAME_
            my @pn = split(/\t/, $mutations[$fn]);
            #_BAD_NAME_
            $t = int( (($fn + $fs)/ 2 ) + 0.5);    
            #_BAD_NAME_
            my @pt = split(/\t/, $mutations[$t]);

            if ($line[1] > $pt[1]) {
                $fs = $t;
            } elsif ($line[2] < $pt[1]) {
                $fn = $t;
            }
            $d+=1; 
        }

        for (my $i = $fs; $i <= $#mutations; $i ++) {
            my @mut  = split(/\t/, $mutations[$i]);
            my $type = assign_type($mut[4], $mut[5]);
            if ($mut[1] >= $line[1] && $mut[1] <= $line[2]) {
                $wgbmr{$type} += 1;
            } elsif ($mut[1] > $line[2]) {
                last; # if sorted array
            } elsif ($mut[1] < $line[1]) {
                next;
            }
        }
    }
    close($inputBED);

    return($wgdist, \%wgbmr);
}

# parse_C3D
# Description:
#   Read C3D output file and store information
# Inputs:
#   c3d_file:           C3D output file path
#   chroms_ref:         pointer to @chroms
#   starts_ref:         pointer to @starts 
#   mutations_ref:      pointer to @mutations
#   window:             window to make bounds on either side of the gene
#   threshold:          correlation threshold
# Outputs:
#   SiMES:              total # of each type of mutation within SRE
#   nbmr:               # of mutations within BMR regions
#   cbmr:               BMR coverage
#   nmut:               # of mutations within test regions
#   nid:                test region coverage
#   mutout_file:        output mutation file being written to
sub parse_C3D {
    my $c3d_file      = shift;
    my $chroms_ref    = shift;
    my @chroms        = @$chroms_ref;
    my $starts_ref    = shift;
    my @starts        = @$starts_ref;
    my $mutations_ref = shift;
    my @mutations     = @$mutations_ref;
    my $window        = shift;
    my $threshold     = shift;

    my %SiMES = ();
    my %nbmr  = ();
    my %cbmr  = ();
    my %nmut  = ();
    my %cmut  = ();
    my %nid   = ();

    my $mutout_file = $options{"Output Directory"} . "/" . join(
        ".",
        basename($c3d_file),
        $window,
        $threshold,
        $suffix_mut
    );

    open(my $inputC3D, "<", $c3d_file) or die "Could not open $c3d_file!\n";
    open(my $mutout, ">", $mutout_file) or die "Could not open $mutout_file!\n";
    #_EXPLANATION_
    while (<$inputC3D>) {
        chomp();
        my @line = split(/\t/);

        if ($line[0] eq "COORD_1") {
            next;
        }
        if ($line[2] eq "NA") { # attempt to fix
            next;
        }

        my @prox          = split(/[:-]/, $line[0]);         # current gene of interest
        my @dist          = split(/[:-]/, $line[1]);         # distal region of comparison
        my $gene          = join(":", $line[0], $line[4]);   # gene region
        my $coverage      = $dist[2] - $dist[1];
        my @window_bounds = ($prox[1] - $window, $prox[2] + $window);  
        
        print $mutout "@dist\t$line[2]\n";

        if (!exists $SiMES{$gene}) {
            $SiMES{$gene} = 0;

            my %mtypes = (
                ACTG  => 0,
                AGTC  => 0,
                ATTA  => 0,
                CAGT  => 0,
                CGGC  => 0,
                CTGA  => 0,
                OTHER => 0,
            );
            my %btypes = (
                ACTG  => 0,
                AGTC  => 0,
                ATTA  => 0,
                CAGT  => 0,
                CGGC  => 0,
                CTGA  => 0,
                OTHER => 0,
            );
            $nbmr{$gene} = \%btypes;
            $nmut{$gene} = \%mtypes; 

            my $id = [];
            $nid{$gene} = $id;
            print $mutout "\n\n---$gene---\n\n";
        }       

        my $start = 0;
        my $fs    = 0;
        my $fn    = 0;

        for (my $i = 0; $i <= $#chroms; $i++) {
            if ($dist[0] eq $chroms[$i]) {
                if ($i != $#chroms) {
                    $fs = $starts[$i];
                    $fn = $starts[$i + 1];
                } else {
                    $fs = $starts[$i];
                    $fn = $#mutations;
                }
            }
        }

        my $d = 0;
        my $t = 0;
        my $s = 10;
        while ($d < $s) {
            my @ps = split(/\t/, $mutations[$fs]);
            my @pn = split(/\t/, $mutations[$fn]);

            $t = int( (($fn + $fs)/ 2) + 0.5);    
            my @pt = split(/\t/, $mutations[$t]);

            if ($dist[1] > $pt[1]) {
                $fs = $t;
            } elsif ($dist[2] < $pt[1]) {
                $fn = $t;
            }
            $d += 1;     
        }

        for (my $i = $fs; $i <= $fn; $i++) {
            my @mut  = split(/\t/,$mutations[$i]);
            my $type = assign_type($mut[4],$mut[5]);

            if ($dist[1] >= $window_bounds[0] && $dist[2] <= $window_bounds[1]) {
                if ($mut[1] >= $dist[1] && $mut[1] <= $dist[2]) {
                    if ($line[2] >= $threshold) {
                        my %lnid = map { $_ => 1 } @{$nid{$gene}};
                        if (!exists $lnid{$mut[3]}) {   
                            print $mutout "@dist\t$line[2]\t$gene\t@mut\tTEST\n";
                            $nmut{$gene}{$type} += 1;
                            push(@{$nid{$gene}}, $mut[3]);
                        } else {
                            print $mutout "@dist\t$line[2]\t$gene\t@mut\tTEST_OMITTED\n";
                        }
                    } else {
                        print $mutout "@dist\t$line[2]\t$gene\t@mut\tlBMR\n";
                        $nbmr{$gene}{$type} += 1;
                    }
                } elsif ($mut[1] > $dist[2]) {
                    last; # if sorted array
                } elsif ($mut[1] < $dist[1]) {
                    next;
                }
            }
        }

        if ($dist[1] >= $window_bounds[0] && $dist[2] <= $window_bounds[1]) {
            if ($line[2] >= $threshold) {
                $SiMES{$gene} += 1;
                $cmut{$gene}  += $coverage;
            } else {
                $cbmr{$gene} += $coverage;
            }
        }
    }
    close($inputC3D);
    close($mutout);

    return(\%SiMES, \%nbmr, \%cbmr, \%nmut, \%cmut, \%nid, $mutout_file);
}

# calculate
# Description:
#   Calculate mutations rates and p-values
# Inputs:
#   output_filename:    MuSE output file name
#   wgdist:             distance of genome analyzed in reference BED (for calculating global mutation rates)
#   wgbmr_ref:          pointer to %wgbmr
#   SiMES_ref:          pointer to %SiMES
#   nbmr_ref:           pointer to %nbmr
#   cbmr_ref:           pointer to %cbmr
#   nmut_ref:           pointer to %nmut
#   cmut_ref:           pointer to %cmut
#   nid_ref:            pointer to %nid
# Outputs:
#   None
sub calculate {
    my $output_filename = shift;
    my $wgdist          = shift;
    my $wgbmr_ref       = shift;
    my %wgbmr           = %$wgbmr_ref;
    my $SiMES_ref       = shift;
    my %SiMES           = %$SiMES_ref;
    my $nbmr_ref        = shift;
    my %nbmr            = %$nbmr_ref;
    my $cbmr_ref        = shift;
    my %cbmr            = %$cbmr_ref;
    my $nmut_ref        = shift;
    my %nmut            = %$nmut_ref;
    my $cmut_ref        = shift;
    my %cmut            = %$cmut_ref;
    my $nid_ref         = shift;
    my %nid             = %$nid_ref;
    my $pm;

    open(my $output, ">", $output_filename) or die "Could not open $output_filename\n";
    print $output $header_line . "\n";

    for my $gene (keys %nmut) {
        my $pid;
        my @pvals  = ();
        my @lpvals = ();
        my @occur  = ();
        my @loccur = ();
        my @rate   = ();
        my @lrate  = ();
        my $bmr;
        my $lbmr;
        print $output "$gene\t";
        
        for my $mtype (sort keys %{$nmut{$gene}}) {
            my $p;
            my $lp;

            if (!defined $cmut{$gene}) {
                $cmut{$gene} = 0;
            }
            if (!defined $cbmr{$gene}) {
                $cbmr{$gene} = 0;
            }

            my $bmr = 0;
            if (0 != $wgdist) {
                $bmr = $wgbmr{$mtype}/$wgdist;
            }
            if ($cbmr{$gene} != 0) {
                $lbmr = $nbmr{$gene}{$mtype}/$cbmr{$gene};
                if ($lbmr == 0) { # get from all DHS sites
                    $lbmr = $bmr;
                }
            } else {
                $lbmr = $bmr;
            }
            
            if ($cmut{$gene} != 0 && $nmut{$gene}{$mtype} > 0) {
                $p  = binomial($nmut{$gene}{$mtype}, $cmut{$gene}, $bmr);
                $lp = binomial($nmut{$gene}{$mtype}, $cmut{$gene}, $lbmr);
            } else {
                $p  = "NA";
                $lp = "NA";
            }
            
            print $output "$mtype\t$bmr\t$lbmr\t$cmut{$gene}\t$nmut{$gene}{$mtype}\t$p\t$lp\t";
            if ($nmut{$gene}{$mtype} > 0) {
                push(@pvals,$p);
                push(@lpvals,$lp);
            }
            push(@occur, $nmut{$gene}{$mtype});
            push(@rate, $bmr);
            push(@lrate, $lbmr);
            
        }

        my $sg = 0;
        my $pcomb;
        my $lpcomb;
        
        if (!@pvals) {
            $pcomb  = "NA";
            $lpcomb = "NA";
        } else {
            $pcomb  = fisher_combine_pval(\@pvals);
            $lpcomb = fisher_combine_pval(\@lpvals);
        }
        print $output "$pcomb\t$lpcomb\t$SiMES{$gene}\n";
    }
    close($output);
}

# split_c3d
# Description:
#   Split C3D output file for cluster processing
# Inputs:
#   c3d_file:       original C3D file to split
# Outputs:
#   file_count:     number of files generated by the split
#   file_list:      pointer to array containing output files
sub split_c3d {
    my $c3d_file = shift;
    my $file_count = 1;
    my $out_filename = $options{"Output Directory"} . "/" . join("_", $prefix_muse, basename($c3d_file), $file_count);
    my $out_current_region = "";
    my $out_filelength = 0;
    my @file_list = ();

    open(my $f_in, "<", $c3d_file) or die "Could not open $c3d_file!\n";
    open(my $f_out, ">", $out_filename) or die "Could not open $out_filename\n";
    while (<$f_in>) {
        my $readline = $_;
        my @splitline = split(/\t/);

        # if onto new block in C3D file, and current file is longer than the max length
        if (($splitline[0] ne $out_current_region) &&
            ($out_filelength > $max_lines_split)) {
            # close file and start new one
            close($f_out);
            push(@file_list, $out_filename);
            push(@tmp_filelist, $out_filename);

            $file_count        += 1;
            $out_filename       = $options{"Output Directory"} . "/" . join("_", $prefix_muse, basename($c3d_file), $file_count);
            $out_current_region = $splitline[0];
            $out_filelength     = 0;
            open($f_out, ">", $out_filename) or die "Could not open $out_filename";
        } elsif ($splitline[0] ne $out_current_region) {
            $out_current_region = $splitline[0];
        }

        print $f_out $readline;
        $out_filelength += 1;
    }
    close($f_in);
    close($f_out);
    push(@file_list, $out_filename);
    push(@tmp_filelist, $out_filename);
    
    return($file_count, \@file_list);
}


### Main ######################################################################
unless (caller) {
    parse_args();
    print("Runtime parameters:\n");
    foreach my $opt_key (sort(keys(%options))) {
        print("\t$opt_key:\t$options{$opt_key}\n");
    }

    # Parse mutations file
    print("Reading mutation file\n");
    my ($mutations_ref, $starts_ref, $chroms_ref) = parse_mutations($options{"Mutations File"});
    print("....Done\n");

    # Parse reference BED file
    print("Reading BED file\n");
    my ($wgdist, $wgbmr_ref) = parse_reference($options{"DHS File"}, $chroms_ref, $starts_ref, $mutations_ref);
    print("....Done\n");

    # Run in parallel (default)
    if ($options{"Parallel"}) {
        my $pm = Parallel::ForkManager->new($max_processes);
        my @muse_filelist = ();
        my @mut_filelist  = ();

        # need to run this so that temp files are properly pushed into arrays (this is just a function of how P:FM works)
        $pm->run_on_finish(sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $pm_filelist_ref) = @_;
            push(@muse_filelist, $pm_filelist_ref->[0]);
            push(@tmp_filelist, $pm_filelist_ref->[0]);
            push(@mut_filelist, $pm_filelist_ref->[1]);
            push(@tmp_filelist, $pm_filelist_ref->[1]);
        });

        # Split C3D file into parts
        print("Splitting C3D file\n");
        my ($file_count, $file_list_ref) = split_c3d($options{"C3D File"});
        print("....Done\n");

        print("Parallel C3D reading and calculating\n");
        foreach my $f (@$file_list_ref) {
            my $pid = $pm->start() and next;

            # Parse C3D file
            my ($SiMES_ref, $nbmr_ref, $cbmr_ref, $nmut_ref, $cmut_ref, $nid_ref, $mut_file) = 
                parse_C3D(
                    $f,
                    $chroms_ref,
                    $starts_ref,
                    $mutations_ref,
                    $options{"Window"},
                    $options{"Threshold"}
                );

            # Run binomial hypothesis test calculations
            my $muse_file = $options{"Output Directory"} . "/" . join(
                ".",
                basename($f),
                $options{"Window"},
                $options{"Threshold"},
                $suffix_muse
            );
            calculate(
                $muse_file,
                $wgdist,
                $wgbmr_ref,
                $SiMES_ref,
                $nbmr_ref,
                $cbmr_ref,
                $nmut_ref,
                $cmut_ref,
                $nid_ref
            );

            my @pm_filelist = ($muse_file, $mut_file);
            $pm->finish(0, \@pm_filelist);
        }
        $pm->wait_all_children();
        print("....Done\n");


        # Cleanup temporary files, collect into single files
        print("Cleaning up\n");
        my $muse_file = $options{"Output Directory"} . "/" . join(
            ".",
            basename($options{"C3D File"}),
            $options{"Window"},
            $options{"Threshold"},
            $suffix_muse
        );
        my $mut_file = $options{"Output Directory"} . "/" . join(
            ".",
            basename($options{"C3D File"}),
            $options{"Window"},
            $options{"Threshold"},
            $suffix_mut
        );

        # Collect MUSE files
        open(my $f_museout, ">", $muse_file) or die "Could not open $muse_file\n";
        print $f_museout $header_line . "\n";
        foreach my $f (@muse_filelist) {
            open(my $f_musein, "<", $f) or die "Could not open $f\n";
            my $readline = <$f_musein>; # skip header line produced by "calculate" function
            while ($readline = <$f_musein>) {
                print $f_museout $readline;
            }
            close($f_musein);
        }
        close($f_museout);

        # Collect MUT files
        open(my $f_mutout, ">", $mut_file) or die "Could not open $mut_file\n";
        foreach my $f (@mut_filelist) {
            open(my $f_mutin, "<", $f) or die "Could not open $f\n";
            while (my $readline = <$f_mutin>) {
                print $f_mutout $readline;
            }
            close($f_mutin);
        }
        close($f_mutout);
        
        remove_tree(@tmp_filelist);
        print("....Done\n");

    # Run in series
    } else {
        # Parse C3D file
        my ($SiMES_ref, $nbmr_ref, $cbmr_ref, $nmut_ref, $cmut_ref, $nid_ref, $mut_file) = 
            parse_C3D(
                $options{"C3D File"},
                $chroms_ref,
                $starts_ref,
                $mutations_ref,
                $options{"Window"},
                $options{"Threshold"}
            );

        # Run binomial hypothesis test calculations
        my $output_file = $options{"Output Directory"} . "/" . join(
            ".",
            basename($options{"C3D File"}),
            $options{"Window"},
            $options{"Threshold"},
            $suffix_muse
        );
        calculate(
            $output_file,
            $wgdist,
            $wgbmr_ref,
            $SiMES_ref,
            $nbmr_ref,
            $cbmr_ref,
            $nmut_ref,
            $cmut_ref,
            $nid_ref
        );

        # Cleanup temporary files
        print("Cleaning up\n");
        remove_tree(@tmp_filelist);
        print("....Done\n");
    }
}


__END__

=head1 NAME

MuSE: calculating significantly mutated regions

=head1 SYNOPSIS

perl MuSE.pm -m <mutations BED> -c <C3D output file> -d <DHS file> -w <window> -t <threshold> [options]

Required Parameters:

    -m,--mutations          BED detail file containing mutations information
    -c,--c3d                Output file from C3D, defining sets of regulatory elements
    -d,-dhs                 DNase Hypersensitivity Sites BED file
    -w,--window             Window (in bp) to define "local" vs "global" mutations. Window must be <= window size used to generate the C3D file (Default: 500kb)
    -t,--threshold          Minimum correlation threshold value (Default: 0.7)

Options:

    -h, --help              Brief help message
    --man                   Man page with full documentation
    -o, --output            Directory path for output files
    --no-parallel           Perform local calculations in series (parallel by default)