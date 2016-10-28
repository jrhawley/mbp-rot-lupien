#!/usr/bin/perl

### Imports ###################################################################
use strict;
use warnings;
use Statistics::R;

# load peak file used for ifc calculate mutation rates...


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

#_COMMENT_: don't make your own function, use prebuilt ones
#_COMMENT_: see http://stackoverflow.com/questions/10701210/how-to-find-maximum-and-minimum-value-in-an-array-of-integers-in-perl for example
#minimum number in an array
sub my_min {
    my $min = shift;
    my $next;
    while(@_){
       $next = shift;
       $min = $next if($next < $min);
    }
    return $min;
}

# assign_type
# Description:
#   Assign a generic type of mutation given the reference and mutant alleles
# Inputs:
#   a: reference allele
#   b: mutant allele
# Outputs:
#   type: type of mutation
sub assign_type{
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
    #_COMMENT_: could also rewrite this function using the %Compl hash

    return($type);
} 

### Main ######################################################################
# Complement for flipping strands
#_NOT_USED_
my %Compl = (
    A => "T",
    C => "G",
    G => "C",
    T => "A"
);


# For calculating the total number of mutations per individual
#_NOT_USED_
my %totalMut_ID=();

# Calculate the total number of each type of mutation
# and total within regulatory elements (IFC).

my %SiMES=();

#_BAD_NAME_
#_NOT_USED_
my %mutations=(); # mutations by chromosome / Type / ID???

# Sorted mutation file
my $inputMUT;
my $inputMUTfile=$ARGV[0];
open($inputMUT, "<", $inputMUTfile) or die "Could not open $inputMUTfile!\n";

#IFC output file
my $inputIFC;
my $inputIFCfile=$ARGV[1];
open($inputIFC, "<", $inputIFCfile) or die "Could not open $inputIFCfile!\n";

#Bed file used for the IFC analysis -- will use genomewide mutation rate if no mutations within selected window
my $inputBED;
my $inputBEDfile=$ARGV[2];
open($inputBED, "<", $inputBEDfile) or die "Could not open $inputBEDfile!\n";

# +/- window size -- Must be equal or smaller than the IFC window size used
my $window=$ARGV[3];

# threshold for SiME / SRE --- ie. elements connected by C3D/IFC will be connect at this threshold
my $thres=$ARGV[4];

my $output;
my $append="MUSE6";
my $outputfile=join(".", $inputIFCfile, $window, $thres, $append);
open($output, ">", $outputfile) or die "Could not open $outputfile\n";

my $mutout;
my $append2="MUT6";
my $mutoutfile=join(".", $inputIFCfile, $window, $thres, $append2);
open($mutout, ">", $mutoutfile) or die "Could not open $mutoutfile\n";

my @mutations=();
my %landmarks=();
my @chroms=();
my @starts=();

my $n=-1; # array position

while(<$inputMUT>){
    chomp();

    $n+=1;
    push(@mutations, $_);

    my @line=split(/\t/);

    my $chr=$line[0];
    my $bp=$line[2];
    my $id=$line[3];
    my $A1=$line[4];
    my $A2=$line[5];
    #_NOT_USED_
    my $mtype;

    #_NOT_USED_
    my $info=join(":", $bp, $A1, $A2, $id);

    #_EXPLANATION_
    if(!exists $landmarks{$chr}){
        $landmarks{$chr}=$n;
        push(@starts, $n);
        push(@chroms, $chr);
    }  

}
close($inputMUT);

#_COMMENT_: add this to the previous while loop, don't go through all the lines again
my %mutP=();
for my $mutation (@mutations){
    my @mut=split(/\t/,$mutation);

    #_EXPLANATION_
    if(!exists $mutP{$mut[3]}){
        $mutP{$mut[3]}=1;
    }else{
        $mutP{$mut[3]}+=1;
    }
}

print("Done InputMUT reading");


my %mgene=();
my %migene=();
my %mexon=();
my $mg=0;
my $mi=0;
my %pmgene=();
my %pmigene=();

my %wgbmr=(
    "ACTG" => 0,
    "AGTC" => 0,
    "ATTA" => 0,
    "CAGT" => 0,
    "CGGC" => 0,
    "CTGA" => 0,
    "OTHER" => 0,
);

my $wgdist=0;

#_EXPLANATION_
while(<$inputBED>){
    chomp();
    my @line=split(/\t/);
    my $dist=$line[2] - $line[1];
    $wgdist+= $dist;

    my $start=0;
    my $fs=0;
    my $fn=0;

    for(my $i=0; $i<=$#chroms; $i++){
        if($line[0] eq $chroms[$i]){
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

        if($line[1] > $pt[1]){
            $fs=$t;
        }elsif($line[2] < $pt[1]){
            $fn=$t;
        }
        $d+=1;     
        #   print "$fs\t$fn\n@ps\n@pn\n";
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

for my $type (keys %wgbmr){
    print "$type $wgbmr{$type} $wgdist\n";
}

my %nbmr=(); #mut within bmr regions
my %cbmr=(); #bmr coverage
my %nmut=(); #mut within test regions
my %cmut=(); #test region coverage
my %tmut=();
my %nid=();

#_EXPLANATION_
while(<$inputIFC>){
    chomp();
    #_COMMENT_: why splitting by space and not \t?
    my @line=split(/ /);

    if($line[0] eq "COORD_1"){
        next;
    }
    if($line[2] eq "NA"){ # attempt to fix
        next;
    }

    my @prox=split(/:/,$line[0]);
    my @dist=split(/:/,$line[1]);

    my $gene=join(":",$line[0],$line[4]); 
    my $cov=$dist[2] - $dist[1];  

    my $lwind=$prox[1] - $window;
    my $uwind=$prox[2] + $window;   
    #print $mutout "@dist\t$line[2]\n";

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
        # print $mutout "\n\n---$gene---\n\n";
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
        #print $mutout "\n---$fs\t$fn\t@ps\t@pn\n"; 
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
    #}

    if($dist[1] >= $lwind && $dist[2] <= $uwind){
        if($line[2] >= $thres){
            $SiMES{$gene}+=1;
            $cmut{$gene}+=$cov;
        }else{
            $cbmr{$gene}+=$cov;
        }
    }
}
close($inputIFC);

#header for output file
#_COMMENT_: use join instead of \t all the time
print $output "GENE\t";
print $output "M_TYPE_1\tM_TYPE_1_BMR\tM_TYPE_1_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_1_N_MUT\tM_TYPE_1_P\tM_TYPE_1_P_LOCAL\t";
print $output "M_TYPE_2\tM_TYPE_2_BMR\tM_TYPE_2_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_2_N_MUT\tM_TYPE_2_P\tM_TYPE_2_P_LOCAL\t";
print $output "M_TYPE_3\tM_TYPE_3_BMR\tM_TYPE_3_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_3_N_MUT\tM_TYPE_3_P\tM_TYPE_3_P_LOCAL\t";
print $output "M_TYPE_4\tM_TYPE_4_BMR\tM_TYPE_4_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_4_N_MUT\tM_TYPE_4_P\tM_TYPE_4_P_LOCAL\t";
print $output "M_TYPE_5\tM_TYPE_5_BMR\tM_TYPE_5_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_5_N_MUT\tM_TYPE_5_P\tM_TYPE_5_P_LOCAL\t";
print $output "M_TYPE_6\tM_TYPE_6_BMR\tM_TYPE_6_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_6_N_MUT\tM_TYPE_6_P\tM_TYPE_6_P_LOCAL\t";
print $output "M_TYPE_7\tM_TYPE_7_BMR\tM_TYPE_7_BMR_LOCAL\tSIME_LENGTH\tM_TYPE_7_N_MUT\tM_TYPE_7_P\tM_TYPE_7_P_LOCAL\t";
print $output "FISHER_P FISHER_P_LOCAL N_DHS\n";

#_EXPLANATION_
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
            $lbmr=$nbmr{$gene}{$mtype}/$cbmr{$gene};
            if($lbmr == 0){ # get from all DHS sites
                $lbmr=$bmr;
            }
        }else{
            $lbmr=$bmr;
        }
        print "\n$nmut{$gene}{$mtype}\t$cmut{$gene}\t$bmr\t$lbmr\n";
        
        if($cmut{$gene} != 0 && $nmut{$gene}{$mtype} > 0){
            $p=binomial($nmut{$gene}{$mtype},$cmut{$gene},$bmr);
            $lp=binomial($nmut{$gene}{$mtype},$cmut{$gene},$lbmr);
        }else{
            $p="NA";
            $lp="NA";
        }
        print "$gene $mtype $cbmr{$gene} $nbmr{$gene}{$mtype} --- $bmr / $lbmr --- $cmut{$gene} $nmut{$gene}{$mtype} $p $SiMES{$gene}\n";
        
        print $output "$mtype\t$bmr\t$lbmr\t$cmut{$gene}\t$nmut{$gene}{$mtype}\t$p\t$lp\t";
        if($nmut{$gene}{$mtype} > 0){
            push(@pvals,$p);
            push(@lpvals,$lp);
        }
        print "PVALS=@pvals\n";
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
    print "\n**** $gene $pcomb $lpcomb\n\n";
    print $output "$pcomb\t$lpcomb\t$SiMES{$gene}\n";
}

close($output);



