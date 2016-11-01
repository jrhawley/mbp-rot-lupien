#!/usr/bin/perl

### Imports ###################################################################
use strict;
use warnings;
use Statistics::R;
use Getopt::Long;
use Pod::Usage;

### Global Variables ##########################################################
my $man              = 0;
my $help             = 0;
my $inputMUTfile     = "";
my $inputC3Dfile     = "";
my $inputBEDfile     = "";
my $window           = 0;
my $thres            = 0;
my $suffix_muse      = "MUSE6";
my $suffix_mut       = "MUT6";
my $output_directory = "";

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
sub binomial{
    my $x=shift;
    my $n=shift;
    my $p=shift;

    # start R instance
    my $R = Statistics::R->new();

    # set variables
    $R->set('x', $x);
    $R->set('n', $n);
    $R->set('p', $p);

    # perform test
    $R->run(q`b <- binom.test(x, n, p, alternative="greater" )`); # may want to make greater than
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
sub fisher_combine_pval{
    #_COMMENT_: why not initially cast as array?
    my $pvals=shift;

    my $R = Statistics::R->new();

    # combine values
    $R->set('pvals', [@{$pvals}]);
    #$R->run(q`chisq <- -2*sum(log(pvals,10))`);
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
    my $a=shift;
    my $b=shift;
    my $type;

    if(($a eq "A" && $b eq "C") || ($a eq "T" && $b eq "G")){
       $type="ACTG";
    }elsif(($a eq "A" && $b eq "G") || ($a eq "T" && $b eq "C")){
       $type="AGTC";
    }elsif(($a eq "A" && $b eq "T") || ($a eq "T" && $b eq "A")){
       $type="ATTA";
    }elsif(($a eq "C" && $b eq "A") || ($a eq "G" && $b eq "T")){
       $type="CAGT";
    }elsif(($a eq "C" && $b eq "G") || ($a eq "G" && $b eq "C")){
       $type="CGGC";
    }elsif(($a eq "C" && $b eq "T") || ($a eq "G" && $b eq "A")){
       $type="CTGA";
    }else{
       $type="OTHER";
    }
    #_COMMENT_: add "NONE" mutation, just for robustness
    return($type);
}

# parse_args
# Description:
#   parse command line input arguments
sub parse_args {
    GetOptions(
        'm|mutations=s' => \$inputMUTfile,
        'c|c3d=s'       => \$inputC3Dfile,
        'b|bed=s'       => \$inputBEDfile,
        'w|window=i'    => \$window,
        't|threshold=f' => \$thres,
        'help|?'        => \$help,
        'man'           => \$man,
        'o|output:s'    => \$output_directory
    ) or pod2usage(2);
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage(-exitval => 0, -verbose => 2);
    }
    pod2usage({ -message => q{Mandatory argument '-m' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $inputMUTfile;
    pod2usage({ -message => q{Mandatory argument '-c' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $inputC3Dfile;
    pod2usage({ -message => q{Mandatory argument '-b' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $inputBEDfile;
    pod2usage({ -message => q{Mandatory argument '-w' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $window;
    pod2usage({ -message => q{Mandatory argument '-t' is missing},
                 -exitval => 1,
                 -verbose => 1
    }) unless $thres;
}

### Main ######################################################################
parse_args();
print("Runtime parameters:\n");
print(join("\n\t", "\t$inputMUTfile", $inputC3Dfile, $inputBEDfile, $window, $thres)."\n");

# Calculate the total number of each type of mutation
# and total within regulatory elements (C3D).
my %SiMES=();

# Sorted mutation file
my $inputMUT;

#C3D output file
my $inputC3D;

#Bed file used for the C3D analysis -- will use genomewide mutation rate if no mutations within selected window
my $inputBED;

# +/- window size -- Must be equal or smaller than the C3D window size used

# threshold for SiME / SRE --- ie. elements connected by C3D will be connect at this threshold

my $output;
my $outputfile = join(".", $inputC3Dfile, $window, $thres, $suffix_muse);

my $mutout;
my $mutoutfile = join(".", $inputC3Dfile, $window, $thres, $suffix_mut);

my @mutations = ();
my %landmarks = ();
my @chroms    = ();
my @starts    = ();

my $n=-1; # array position

print("Reading mutation file\n");
open($inputMUT, "<", $inputMUTfile) or die "Could not open $inputMUTfile!\n";
while (<$inputMUT>) {
    chomp();

    $n+=1;
    push(@mutations, $_);

    # extract information from mutations file
    my @line = split(/\t/);
    my $chr  = $line[0];
    my $bp   = $line[2];
    my $id   = $line[3];
    my $A1   = $line[4];
    my $A2   = $line[5];

    #_EXPLANATION_
    if (!exists $landmarks{$chr}) {
        $landmarks{$chr} = $n;
        push(@starts, $n);
        push(@chroms, $chr);
    }  

}
close($inputMUT);

#_COMMENT_: add this to the previous while loop, don't go through all the lines again
my %mutP=();
for my $mutation (@mutations) {
    my @mut = split(/\t/, $mutation);

    #_EXPLANATION_
    if (!exists $mutP{$mut[3]}) {
        $mutP{$mut[3]} = 1;
    } else {
        $mutP{$mut[3]} += 1;
    }
}
print("Finished reading\n");


my %mgene=();
my %migene=();
my %mexon=();
my $mg=0;
my $mi=0;
my %pmgene=();
my %pmigene=();
my $wgdist=0;
my %wgbmr=(
    "ACTG" => 0,
    "AGTC" => 0,
    "ATTA" => 0,
    "CAGT" => 0,
    "CGGC" => 0,
    "CTGA" => 0,
    "OTHER" => 0,
);

print("Reading BED file\n");
open($inputBED, "<", $inputBEDfile) or die "Could not open $inputBEDfile!\n";
#_EXPLANATION_
while (<$inputBED>) {
    chomp();
    my @line = split(/\t/);
    my $dist = $line[2] - $line[1];
    $wgdist+= $dist;

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

    for (my $i=$fs; $i <= $#mutations; $i ++){
        my @mut=split(/\t/,$mutations[$i]);
        my $type=assign_type($mut[4],$mut[5]);
        if($mut[1] >= $line[1] && $mut[1] <= $line[2]){
            $wgbmr{$type}+=1;
        }elsif( $mut[1] > $line[2] ){
            last; # if sorted array
        }elsif( $mut[1] < $line[1] ){
            next;
        }
    }
}
close($inputBED);
print("Finished reading\n");

my %nbmr=(); #mut within bmr regions
my %cbmr=(); #bmr coverage
my %nmut=(); #mut within test regions
my %cmut=(); #test region coverage
my %tmut=();
my %nid=();

print("Reading C3D output\n");
open($inputC3D, "<", $inputC3Dfile) or die "Could not open $inputC3Dfile!\n";
open($mutout, ">", $mutoutfile) or die "Could not open $mutoutfile\n";
#_EXPLANATION_
while (<$inputC3D>) {
    chomp();
    my @line=split(/\t/);

    if($line[0] eq "COORD_1"){
        next;
    }
    if($line[2] eq "NA"){ # attempt to fix
        next;
    }

    my @prox=split(/[:-]/,$line[0]);
    my @dist=split(/[:-]/,$line[1]);

    my $gene=join(":",$line[0],$line[4]);

    my $cov=$dist[2] - $dist[1];  

    my $lwind=$prox[1] - $window;
    my $uwind=$prox[2] + $window;   
    print $mutout "@dist\t$line[2]\n";

    if(!exists $SiMES{$gene}){
        $SiMES{$gene}=0;

        my %mtypes=(
            ACTG => 0,
            AGTC => 0,
            ATTA => 0,
            CAGT => 0,
            CGGC => 0,
            CTGA => 0,
            OTHER => 0,
        );
        my %btypes=(
            ACTG => 0,
            AGTC => 0,
            ATTA => 0,
            CAGT => 0,
            CGGC => 0,
            CTGA => 0,
            OTHER => 0,
        );
        my $btype=\%btypes;
        my $mtype=\%mtypes;
        $nbmr{$gene}=$btype;
        $nmut{$gene}=$mtype; 

        my $id=[];
        $nid{$gene}=$id;
        print $mutout "\n\n---$gene---\n\n";
    }       

    my $start=0;
    my $fs=0;
    my $fn=0;

    for(my $i=0; $i<=$#chroms; $i++){
        if($dist[0] eq $chroms[$i]){
            if($i != $#chroms){
                $fs=$starts[$i];
                $fn=$starts[$i + 1];
            }else{
                $fs=$starts[$i];
                $fn=$#mutations;
            }
        }
    }

    my $d = 0;
    my $t = 0;
    my $s = 10;
    while( $d < $s ){
        my @ps=split(/\t/,$mutations[$fs]);
        my @pn=split(/\t/,$mutations[$fn]);

        $t=int( (($fn + $fs)/ 2 )+0.5);    
        my @pt=split(/\t/,$mutations[$t]);

        if($dist[1] > $pt[1]){
            $fs=$t;
        }elsif($dist[2] < $pt[1]){
            $fn=$t;
        }
        $d+=1;     
        print $mutout "\n---$fs\t$fn\t@ps\t@pn\n"; 
    }

    for (my $i=$fs; $i <= $fn; $i++){
        my @mut=split(/\t/,$mutations[$i]);
        my $type=assign_type($mut[4],$mut[5]);
        if($dist[1] >= $lwind && $dist[2] <= $uwind){
                if($mut[1] >= $dist[1] && $mut[1] <= $dist[2]){
                if($line[2] >= $thres){
                    #push(@{$nid{$gene}}, $mut[3]);
                    my %lnid = map { $_ => 1 } @{$nid{$gene}};
                    if(!exists $lnid{$mut[3]}){   
                        print $mutout "@dist\t$line[2]\t$gene\t@mut\tTEST\n";
                        $nmut{$gene}{$type}+=1;
                        push(@{$nid{$gene}}, $mut[3]);
                    }else{
                        print $mutout "@dist\t$line[2]\t$gene\t@mut\tTEST_OMITTED\n";
                    }
                }else{ 
                    print $mutout "@dist\t$line[2]\t$gene\t@mut\tlBMR\n";
                    $nbmr{$gene}{$type}+=1;
                }
            }elsif( $mut[1] > $dist[2] ){
                last; # if sorted array
            }elsif( $mut[1] < $dist[1] ){
                next;
            }
        }
    }

    if($dist[1] >= $lwind && $dist[2] <= $uwind){
        if($line[2] >= $thres){
            $SiMES{$gene}+=1;
            $cmut{$gene}+=$cov;
        }else{
            $cbmr{$gene}+=$cov;
        }
    }
}
close($inputC3D);
close($mutout);
print("Finished reading\n");

#header for output file
#_COMMENT_: use join instead of \t all the time
open($output, ">", $outputfile) or die "Could not open $outputfile\n";
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
    "M_TYPE_2",
    "M_TYPE_2_BMR",
    "M_TYPE_2_BMR_LOCAL",
    "SIME_LENGTH",
    "M_TYPE_2_N_MUT",
    "M_TYPE_2_P",
    "M_TYPE_2_P_LOCAL",
    "M_TYPE_2",
    "M_TYPE_2_BMR",
    "M_TYPE_2_BMR_LOCAL",
    "SIME_LENGTH",
    "M_TYPE_2_N_MUT",
    "M_TYPE_2_P",
    "M_TYPE_2_P_LOCAL",
    "M_TYPE_2",
    "M_TYPE_2_BMR",
    "M_TYPE_2_BMR_LOCAL",
    "SIME_LENGTH",
    "M_TYPE_2_N_MUT",
    "M_TYPE_2_P",
    "M_TYPE_2_P_LOCAL",
    "M_TYPE_2",
    "M_TYPE_2_BMR",
    "M_TYPE_2_BMR_LOCAL",
    "SIME_LENGTH",
    "M_TYPE_2_N_MUT",
    "M_TYPE_2_P",
    "M_TYPE_2_P_LOCAL",
    "M_TYPE_2",
    "M_TYPE_2_BMR",
    "M_TYPE_2_BMR_LOCAL",
    "SIME_LENGTH",
    "M_TYPE_2_N_MUT",
    "M_TYPE_2_P",
    "M_TYPE_2_P_LOCAL",
    "FISHER_P",
    "FISHER_P_LOCAL",
    "N_DHS"
);
print $output $header_line . "\n";

#_EXPLANATION_
print("Starting calculations\n");
for my $gene (keys %nmut){
    my @pvals=();
    my @lpvals=();
    my @occur=();
    my @loccur=();
    my @rate=();
    my @lrate=();
    my $bmr;
    my $lbmr;
    print $output "$gene\t";
    
    for my $mtype ( sort keys %{$nmut{$gene}}){
        #my $tmut=$nbmr{$gene}{$mtype} + $nmut{$gene}{$mtype};
        #my $tcov=$cbmr{$gene} + $cmut{$gene};
        my $p;
        my $lp;
        if(!defined $cmut{$gene}){
            $cmut{$gene}=0;
        }
        if(!defined $cbmr{$gene}){
            $cbmr{$gene}=0;
        }
        my $bmr=$wgbmr{$mtype}/$wgdist;
        if($cbmr{$gene} != 0){
            #_COMMENT_: possible illegal division by 0 here
            $lbmr=$nbmr{$gene}{$mtype}/$cbmr{$gene};
            if($lbmr == 0){ # get from all DHS sites
                $lbmr=$bmr;
            }
        }else{
            $lbmr=$bmr;
        }
        # print "\n$nmut{$gene}{$mtype}\t$cmut{$gene}\t$bmr\t$lbmr\n";
        
        if($cmut{$gene} != 0 && $nmut{$gene}{$mtype} > 0){
            $p=binomial($nmut{$gene}{$mtype},$cmut{$gene},$bmr);
            $lp=binomial($nmut{$gene}{$mtype},$cmut{$gene},$lbmr);
        }else{
            $p="NA";
            $lp="NA";
        }
        # print "$gene $mtype $cbmr{$gene} $nbmr{$gene}{$mtype} --- $bmr / $lbmr --- $cmut{$gene} $nmut{$gene}{$mtype} $p $SiMES{$gene}\n";
        
        print $output "$mtype\t$bmr\t$lbmr\t$cmut{$gene}\t$nmut{$gene}{$mtype}\t$p\t$lp\t";
        if($nmut{$gene}{$mtype} > 0){
            push(@pvals,$p);
            push(@lpvals,$lp);
        }
        # print "PVALS=@pvals\n";
        push(@occur,$nmut{$gene}{$mtype});
        push(@rate,$bmr);
        push(@lrate,$lbmr);
    }

    my $sg=0;
    my $pcomb;
    my $lpcomb;
    
    if(!@pvals){
        $pcomb="NA";
        $lpcomb="NA";
    }else{
        $pcomb=fisher_combine_pval(\@pvals);
        $lpcomb=fisher_combine_pval(\@lpvals);
    }
    # print "\n**** $gene $pcomb $lpcomb\n\n";
    print $output "$pcomb\t$lpcomb\t$SiMES{$gene}\n";
}
close($output);
print("Finished calculations\n");

__END__

=head1 NAME

MuSE: calculating significantly mutated regions

=head1 SYNOPSIS

perl MuSE.pl -m <mutations BED> -c <C3D output file> -b <DHS reference BED> -w <window> -t <threshold>

Options:
    --help              Brief help message
    --man               Man page with full documentation
    -o, --output        Directory path for output files