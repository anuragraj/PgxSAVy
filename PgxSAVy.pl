#!/usr/bin/perl

###########################################################################################
# To calculate VAS score for variants identified
# Input: Result files and Spectra file
# OUTPUT : variant Score file
# Created 01 Dec 2022, Anurag
# Updated 12 June 2023, Anurag
###########################################################################################

use strict;
use warnings;

use FindBin qw( $RealBin );
use lib $RealBin;

use Term::ANSIColor;
my $VersionInfo = "1.0"; #version

#local modules
use MODULES::MSFileParserAR;
use MODULES::Scoring;
use MODULES::Isobaric;
use MODULES::MSGF_modif;
use MODULES::BioAnnotation;

# CPAN modules
use InSilicoSpectro::InSilico::MassCalculator;
use InSilicoSpectro::InSilico::Peptide;
InSilicoSpectro::init('InSilicoSpectro/config/insilicodef.xml');

use Statistics::Zed;
use Getopt::Std;

use constant H => 1.007825;

print color('bold blue');
print "\t","="x60,"\n";
print color('bold yellow');
print "\t\t".' ____             ____    ___     __     '."\n";
print "\t\t".'|  _ \ __ ___  __/ ___|  / \ \   / /   _ '."\n";
print "\t\t".'| |_) / _` \ \/ /\___ \ / _ \ \ / / | | |'."\n";
print "\t\t".'|  __/ (_| |>  <  ___) / ___ \ V /| |_| |'."\n";
print "\t\t".'|_|   \__, /_/\_\|____/_/   \_\_/  \__, |'."\n";
print "\t\t".'      |___/                        |___/ '."\n";

print color('bold blue');
print "\t\t\t   Variant Peptides Scoring\n\t  (Annotating the variant peptides using VAS scoring)\n";
print "\t\t\t\tVersion: $VersionInfo\n";

print "\t","="x60,"\n";
print color('reset');

$|++; # autoflush

my $USAGE=<<"EOU";
-------------------------------------------------------------------------------------
Program:    $0 ($VersionInfo)
Purpose:    Annotating the variant peptides using VAS scoring.
            Specify a path for file containing variant peptides details and spectra file.
            
Usage:      $0 [options] <variant_search_file> <MGF_file>

Options:
            -s <int>        : Search result type (default: 1 [Tab-Seperated])
            -e <int>        : Search engine count (default: 1)
            -t <float>      : Fragment tolerance in Da (default: 0.5)
            -f <int>        : Fixed modification (default: 3)
            -v <int>        : Variable modification (default: 1)
            -i <int>        : Instrument type (default: 6)
            -p <int>        : Peak threshold (default: 0)
            -m <int>        : Minimum peak (default: 7)
            -b <True|False> : BioAnnotation check (default: TRUE)
            -a <Filename>   : BioAnnotation database file name (default: homo_sapiens_variation.txt)
            -d <True|False> : Keep duplicates (default: FALSE)
            -o <Filename>   : Output file name (default: InputFile_VAS_score_jobtime.tsv)
            -h              : print this usage statement
-------------------------------------------------------------------------------------
EOU

#command line flags/options
my %options=();
getopts('s:e:t:f:v:i:p:m:b:a:d:o:h', \%options);

#print help if requested
die $USAGE if $options{'h'};

#check for required input
my $InputFile = shift or die "Error: Please provide an input variant search file.\n".$USAGE;
my $spectraFile = shift or die "Error: Please provide an input spectra file.\n".$USAGE;

if (! -e $InputFile) { die  "Error: could not find input file $InputFile!\n".$USAGE; }
if (! -e $spectraFile) { die  "Error: could not find spectra file $spectraFile!\n".$USAGE; }

#Check if paramerters are provided, if not use default values
my $searchresulttype = $options{'s'} || 1;
my $total_se_count = $options{'e'} || 1;
my $frag_tol = $options{'t'} || 0.5;
my $fixed_mod = $options{'f'} || 3;
my $var_mod = $options{'v'} || 1;
my $ins_type = $options{'i'} || 6;
my $peak_threshold = $options{'p'} || 0;
my $minimum_peak = $options{'m'} || 7;
my $keep_duplicates = uc($options{'d'} || 'FALSE');
my $OutputFile = $options{'o'} || '';
my $unit='Da';
my $BioAnno_check = uc($options{'b'} || 'TRUE');
my $BioAnno_file = $options{'a'} || 'homo_sapiens_variation.txt';

print "\nParameters used:\n";
print "Search result type: $searchresulttype\n";
print "Total search engine count: $total_se_count\n";
print "Fragment tolerance: $frag_tol $unit\n";
print "Fixed modification: $fixed_mod\n";
print "Variable modification: $var_mod\n";
print "Instrument type: $ins_type\n";
print "Peak threshold: $peak_threshold\n";
print "Minimum peak: $minimum_peak\n";
print "BioAnnotation: $BioAnno_check\n";
print "Keep duplicates: $keep_duplicates\n";

print "\n";

#check required files
my $required_files1 = "InSilicoSpectro.pm";
my $required_files2 = "homo_sapiens_variation_bioAnnoDB.tsv";
my $required_dir1 = "MODULES";
my $required_dir2 = "InSilicoSpectro";

if (! -e $required_files1) { die  "Error: InSilicoSpectro.pm file missing. Download PgxSAVy again!\n"; }
if (! -d $required_dir1) { die  "Error: MODULES directory missing. Download PgxSAVy again!\n"; }
if (! -d $required_dir2) { die  "Error: InSilicoSpectro directory missing. Download PgxSAVy again!\n"; }

if($BioAnno_check eq "TRUE"){
    if (! -e $BioAnno_file) { die  "Error: could not find BioAnnotation database file \($BioAnno_file\) when BioAnnotation check is true.\nYou can download BioAnnotation database file fom UniProt. Check Documentation!\n"; }
}

#time estimation
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my @abbr = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
$year += 1900;
my $jobtime = "$mday$abbr[$mon]$year"."_".sprintf "%02u%02u%02u", $hour, $min, $sec;
print "Variant Scoring started at ".localtime().".\n";

my $t0=time();

if($OutputFile eq ""){
    $OutputFile = $InputFile;
    $OutputFile =~ s/\.\w+$//;
    $OutputFile .= "_VAS_score_".$jobtime.".tsv";
}

#defining variables
my %scans; #scan list to fetch spectra data 
my %psm_counts; #PSM counts per peptide
#array to store and process input data
my @data_array;

#getting the scan lists
print "Reading the scan list from $InputFile...\n";

#reading the input file
open(IN, $InputFile) or die "Error in reading file $InputFile: $!";
my $header = <IN>;
chomp $header;
$header =~ s/\r//;

open (LOG, ">>PgxSAVyLog.txt") or die "Error in writing file PgxSAVyLog.txt: $!";
print LOG "Job initiated for $InputFile at $jobtime\n";
print LOG "Parameters used: Fragment tolerance: $frag_tol$unit, Fixed modification: $fixed_mod, Variable modification: $var_mod, Instrument type: $ins_type, Peak threshold: $peak_threshold, Minimum peak: $minimum_peak\n";

if($searchresulttype == 3){
    #EuGenoSuite result file
    while (my $line=<IN>)
    {
        chomp $line;
        $line =~ s/\r//;

        #splitting the line to parse the data
        my @arr=split(/\t/,$line);
        
        #defining and assigning variables
        my ($scan,$se,$peptide,$protein_id,$modif_string,$charge_state) = @arr[0,1,3,4,10,6]; 

        my $peptide_length = length($peptide);
        
        #storing the scan list
        $scans{$scan}=undef;

        #counting the number of PSMs per scan
        $psm_counts{$peptide}{$modif_string}++;

        #search engines identification counts
        my $se_counts = 0;
        if($se == 1){
            $se_counts = 3;
        }elsif($se == 2 or $se == 3 or $se == 4){
            $se_counts = 2;
        }elsif($se == 5 or $se == 6 or $se == 7){
            $se_counts = 1;
        }else{
            #any case where the search engine is not identified
            $se_counts = 1;
        }

        # when protein is contaminant or decoy
        if($protein_id =~ /Contaminant/ or $protein_id =~ /DECOY/)
        {
            #skip as it is a contaminant or decoy
        }
        #when protien is variant
        elsif($protein_id =~ /nxp/)
        {
            my @proteins = split(/;nxp/,$protein_id);
            my $wt_peptide = $peptide;
            
            foreach my $protein_id(@proteins){

                my $var_aa = "";
                my $wt_aa = "";
                
                my $var_pos;
                my $var_pos1 = 0;
                my $var_pos2 = 0;
                my $pos2_flag = 0;

                my $protein_pos;

                #parsing variant position from protein ID
                #nxp:NX_Q9H330-4(TP:311-337)(Var:P|326|A;S|332|L)
                $protein_id =~ m/TP:([0-9]+)-([0-9]+)/;
                my ($start,$end) = ($1,$2);
                my ($wt1,$pos1,$mut1,$wt2,$pos2,$mut2);
                if($protein_id =~ m/;/ && $protein_id !~ m/;$/){
                    $protein_id =~ m/Var:([A-Z])\|([0-9]+)\|([A-Z]);([A-Z])\|([0-9]+)\|([A-Z])/;
                    ($wt1,$pos1,$mut1,$wt2,$pos2,$mut2) = ($1,$2,$3,$4,$5,$6);
                    
                    $pos2_flag = 1;
                    $var_aa = $mut1.";".$mut2;
                    $wt_aa = $wt1.";".$wt2;

                    $var_pos1 = ($pos1 - $start)+1;
                    $var_pos2 = ($pos2 - $start)+1;
                    $var_pos = $var_pos1.";".$var_pos2;

                    $protein_pos = $pos1.";".$pos2;
                }
                else
                {
                    $protein_id =~ m/Var:([A-Z])\|([0-9]+)\|([A-Z])/;
                    ($wt1,$pos1,$mut1) = ($1,$2,$3);
                    
                    $var_aa = $mut1;
                    $wt_aa = $wt1;

                    $var_pos1 = ($pos1 - $start)+1;
                    $var_pos = $var_pos1;

                    $protein_pos = $pos1;
                }

                #conditions to check for variant position
                if(not defined $pos1 or not defined $start){
                    print LOG "Warning: variant position not found for $protein_id with $scan.\n";
                    next;
                }
                if($pos1 < $start or $pos1 > $end){
                    print LOG "Warning: variant position mismatch for $protein_id with $scan.\n";
                    next;
                }
                #for those where one index is shifted
                if(($end-$start) == $peptide_length){
                    $start += 1;
                }
            
                #print "VAR1: $var_pos1\tVAR2: $var_pos2\t$peptide_length\t$peptide\t$protein_id\n"; <STDIN>;
                #creating the wild type peptide by substituting the variant position with wild type amino acid
                if($peptide_length >= abs($var_pos1-1) && $peptide_length >= abs($var_pos2-1) && $wt_peptide ne ""){
                    substr($wt_peptide, $var_pos1-1, 1) = $wt1;
                    if($pos2_flag == 1){
                        substr($wt_peptide, $var_pos2-1, 1) = $wt2;
                    }
                }else{
                    print LOG "Warning: variant position mismatch for $protein_id with $scan.\n";
                    next;
                } 

                #defining variant type
                my $var_type = "";
                if($pos2_flag == 1){
                    $var_type = "Double";
                }else{
                    $var_type = "Single";
                }
                
                #create array of arrays for each scan
                my @temp = ($scan,$peptide,$wt_peptide,$var_type,$var_aa,$wt_aa,$var_pos,$protein_pos,$peptide_length,$protein_id,$modif_string,$charge_state,$se_counts);
                push (@data_array, \@temp);

            }
        }else{
            #unknown protein type
        }
    }
}
elsif($searchresulttype == 2){
    #MSGF result output
    while (my $line=<IN>)
    {
        chomp $line;
        $line =~ s/\r//;

        #splitting the line to parse the data
        my @arr=split(/\t/,$line);
        
        #defining and assigning variables #TODO
        my ($scan,$peptide,$protein_id,$modif,$charge_state) = @arr[0,4,6,5,9]; 

        if($protein_id =~ /Contaminant/ or $protein_id =~ /DECOY/)
        {
            #skip as it is a contaminant or decoy
            next;
        }

        my $se_counts = 1;
        my $peptide_length = length($peptide);

        #storing the scan list
        $scans{$arr[0]}=undef;

        #create a modification string
        my $fixed_mod_str = "";
        my $var_mod_str = "";
        my $modif_string = MSGF_modif($peptide,$modif,$fixed_mod,$var_mod);

        #counting the number of PSMs per scan
        $psm_counts{$peptide}{$modif_string}++;

        if($protein_id =~ /nxp/)
        {
            my @proteins = split(/;nxp/,$protein_id);
            my $wt_peptide = $peptide;
            
            foreach my $protein_id(@proteins){

                my $var_aa = "";
                my $wt_aa = "";
                
                my $var_pos;
                my $protein_pos;
                my $var_pos1 = 0;
                my $var_pos2 = 0;
                my $pos2_flag = 0;
                #parsing variant position from protein ID
                #nxp:NX_Q9H330-4(TP:311-337)(Var:P|326|A;S|332|L)
                $protein_id =~ m/TP:([0-9]+)-([0-9]+)/;
                my ($start,$end) = ($1,$2);
                my ($wt1,$pos1,$mut1,$wt2,$pos2,$mut2);
                if($protein_id =~ m/;/ && $protein_id !~ m/;$/){
                    $protein_id =~ m/Var:([A-Z])\|([0-9]+)\|([A-Z]);([A-Z])\|([0-9]+)\|([A-Z])/;
                    ($wt1,$pos1,$mut1,$wt2,$pos2,$mut2) = ($1,$2,$3,$4,$5,$6);
                    
                    $pos2_flag = 1;
                    $var_aa = $mut1.";".$mut2;
                    $wt_aa = $wt1.";".$wt2;

                    $var_pos1 = ($pos1 - $start)+1;
                    $var_pos2 = ($pos2 - $start)+1;
                    $var_pos = $var_pos1.";".$var_pos2;

                    $protein_pos = $pos1.";".$pos2;
                }
                else
                {
                    $protein_id =~ m/Var:([A-Z])\|([0-9]+)\|([A-Z])/;
                    ($wt1,$pos1,$mut1) = ($1,$2,$3);
                    
                    $var_aa = $mut1;
                    $wt_aa = $wt1;

                    $var_pos1 = ($pos1 - $start)+1;
                    $var_pos = $var_pos1;

                    $protein_pos = $pos1;
                }

                #conditions to check for variant position
                if(not defined $pos1 or not defined $start){
                    print LOG "Warning: variant position not found for $protein_id with $scan.\n";
                    next;
                }
                if($pos1 < $start or $pos1 > $end){
                    print LOG "Warning: variant position mismatch for $protein_id with $scan.\n";
                    next;
                }
                #for those where one index is shifted
                if(($end-$start) == $peptide_length){
                    $start += 1;
                }

                #print wt2 if uninitialized #DEBUG
                if($pos2_flag ==1 && not defined $wt2){
                    print "WT Error:$peptide\t$protein_id\n";
                }
            
                #print "VAR1: $var_pos1\tVAR2: $var_pos2\tPL:$peptide_length\tPEP:$peptide\tWT:$wt_peptide\tProt:$protein_id\n"; <STDIN>;
                #creating the wild type peptide by substituting the variant position with wild type amino acid
                if($peptide_length >= abs($var_pos1-1) && $peptide_length >= abs($var_pos2-1) && $wt_peptide ne ""){
                    substr($wt_peptide, $var_pos1-1, 1) = $wt1;
                    if($pos2_flag == 1){
                        substr($wt_peptide, $var_pos2-1, 1) = $wt2;
                    }
                }else{
                    print LOG "Warning: variant position mismatch for $protein_id with $scan.\n";
                    next;
                } 

                #defining variant type
                my $var_type = "";
                if($pos2_flag == 1){
                    $var_type = "Double";
                }else{
                    $var_type = "Single";
                }

                #create array of arrays for each scan
                my @temp = ($scan,$peptide,$wt_peptide,$var_type,$var_aa,$wt_aa,$var_pos,$protein_pos,$peptide_length,$protein_id,$modif_string,$charge_state,$se_counts);
                push (@data_array, \@temp);
            }
        }else{
            #unknown protein type
        }

    }

}elsif($searchresulttype == 1){
    #text result output
    while (my $line=<IN>)
    {
        chomp $line;
        $line =~ s/\r//;

        #splitting the line to parse the data
        my @arr=split(/\t/,$line);

        if(scalar(@arr) != 12){
            print "Error: Incorrect number of columns in the input file.\n";
            print "Expected 12 columns, but found ".scalar(@arr)." columns.\n";
            print "Please check the input file.\n";
            exit;
        }
        
        #defining and assigning variables
        my ($scan,$peptide,$wt_peptide,$var_type,$var_aa,$wt_aa,$var_pos,$protein_pos,$protein_id,$modif_string,$charge_state,$se_counts) = @arr; 

        my $peptide_length = length($peptide);

        #storing the scan list
        $scans{$arr[0]}=undef;

        #counting the number of PSMs per scan
        $psm_counts{$peptide}{$modif_string}++;

        #create array of arrays for each scan
        my @temp = ($scan,$peptide,$wt_peptide,$var_type,$var_aa,$wt_aa,$var_pos,$protein_pos,$peptide_length,$protein_id,$modif_string,$charge_state,$se_counts);
        push (@data_array, \@temp);
    }

}
close IN;

print "Reading the scan details from spectra file $spectraFile...\n";
## Reading spectra from MGF file using spectra_parser
my $ExpSpecRef = spectra_parser( $peak_threshold, $minimum_peak, $spectraFile,\%scans);
my %ExpSpectra = %$ExpSpecRef;
print "Done reading spectra file. ";
print "Total spectra stored = ", scalar keys %ExpSpectra, "\n\n";

#############################################################
#reading results input file
#############################################################
open(IN, $InputFile) or die $!;

print "Calculating scores for spectrum from file $InputFile...\n";

#array to store results
my @results;
my $k=0;
my $inputs_cols = 0;
my $i=0; #for diff spectra scans
my $j=0; # spectrum 2-D for each scan

#while (my $line=<IN>)

#read @data_array and calculate scores
for my $row (0..$#data_array) 
{
    
    my ($scan,$peptide,$wt_peptide,$var_type,$var_aa,$wt_aa,$var_pos,$protein_pos,$peptide_length,$protein_id,$modif_string,$charge_state,$se_counts) = @{ $data_array[$row] };

    my $sum_intensity = 0;
    #hash for the variant position ions
    my %variant_ions = (); #NOT BEING USED CURRENTLY #to be removed from mainscoringcaller

    if (exists $ExpSpectra{$scan}) 
    {
        my @ExpSpec = @{$ExpSpectra{$scan}};
        my $mod_key = $scan;
        my $spec_string=$ExpSpec[5];
        if (not defined $spec_string)
        {
            $k++;
            next;
        }
        my @spec_lines=split("\n",$spec_string);
        pop @spec_lines;
        my @exp=(); #2D array to store spectra data
        my $line_counter=0;
        foreach my $line (@spec_lines)
        {
            my @line_values=split(' ',$line);
            $exp[$line_counter][0]=$line_values[0];
            $exp[$line_counter][1]=$line_values[1];
            $line_counter++;
            
        }
        @spec_lines=();
        for ( my $j = 0 ; $j <= $#exp ; $j++ )
        {
            $sum_intensity += $exp[$j][1];
        }

        #PSM Count
        my $psm_count = $psm_counts{$peptide}{$modif_string};

        #defining variables for global access
        my ($masswiz_score,$fractional_intensity_covered,$rmsd, $bycontinuity_count,$norm_var_score,$avg_var_score,$norm_var_intensity,$avg_var_intensity,$deltaMassWizWT) ;
        my ($wt_masswiz_score,$wt_normMassWizScore,$wt_fractional_intensity_covered,$wt_rmsd, $wt_bycontinuity_count,$wt_norm_var_score,$wt_avg_var_score,$wt_norm_var_intensity,$wt_avg_var_intensity);
        my @decoy_scores;
        my ($second_rank_decoy_score, $normSecond_rank_decoy_score) = (0,0);

        #---------------------------------------------------------------
        ## Calculate the variant score using scoring module
        #---------------------------------------------------------------
        ($masswiz_score,$fractional_intensity_covered,$rmsd, $bycontinuity_count,$norm_var_score,$avg_var_score,$norm_var_intensity,$avg_var_intensity) = MainScoringCaller(\@exp, $peptide, $sum_intensity, $frag_tol, $modif_string, $charge_state, $unit, $ins_type, %variant_ions);
                                
        ##---------------------------------------------------------------
        ## Generate decoy peptides by moving varaint aa to each position
        ## Calculate decoy scores and compare with variant score
        ##----------------------------------------------------------------
        my ($current_var_pos1,$current_var_pos2) = (0,0);
        if($var_type eq "Single"){
            $current_var_pos1 = $var_pos-1;
        }
        else{
            my @var_pos = split(/\;/,$var_pos);
            $current_var_pos1 = $var_pos[0]-1;
            $current_var_pos2 = $var_pos[1]-1;
        }

        if($current_var_pos1 >= $peptide_length or $current_var_pos2 >= $peptide_length){
            #print "Variant position is greater than peptide length for peptide $peptide\n";
            next;
        }

        my $var_pos1 = $current_var_pos1; #for second varinat decoy score calculation

        #loop for each aa position
        for (my $i = 0; $i<$peptide_length; $i++){
            my $vardecoy_peptide = $peptide;
            my $vardecoy_modif_string = $modif_string;
            #move the variant aa to each position
            if($i != $current_var_pos1 and $peptide_length >= abs($current_var_pos1)){      
                my @vardecoy_peptides = split(//,$vardecoy_peptide);    
                ($vardecoy_peptides[$current_var_pos1],$vardecoy_peptides[$i]) = ($vardecoy_peptides[$i],$vardecoy_peptides[$current_var_pos1]);
                if($vardecoy_peptides[$i] eq ""){
                    print "vardecoy_peptide = $scan\t$peptide\t@vardecoy_peptides\t$current_var_pos1\t$i\n";
                    print "AA:$vardecoy_peptides[$current_var_pos1]\t$vardecoy_peptides[$i]\n";
                } 
                $vardecoy_peptide = join('',@vardecoy_peptides);
                my @modrev = split /:/,$vardecoy_modif_string,-1;
                ($modrev[$i+1],$modrev[$current_var_pos1+1]) = ($modrev[$current_var_pos1+1],$modrev[$i+1]);
                $vardecoy_modif_string = join ':',@modrev;
            }
            else{
                next;
            } 

            #skip if the decoy peptide is same as variant peptide
            if ($peptide eq $vardecoy_peptide) {
                next;
            }

            #calculate decoy scores
            my ($dec_masswiz_score,$dec_fractional_intensity_covered,$dec_rmsd, $dec_bycontinuity_count,$norm_var_score,$avg_var_score,$norm_var_intensity,$avg_var_intensity) = MainScoringCaller(\@exp, $vardecoy_peptide, $sum_intensity, $frag_tol, $modif_string, $charge_state, $unit, $ins_type, %variant_ions);

            # Storing decoy scores
            push(@decoy_scores,$dec_masswiz_score);

        }

        #for second variant if present
        if($var_type eq "Double"){
            
            #loop for each aa position
            for (my $i = 0; $i<$peptide_length; $i++){
                my $vardecoy_peptide = $peptide;
                my $vardecoy_modif_string = $modif_string;
                #move the variant aa to each position
                if($i != $current_var_pos2 and $i != ($var_pos1-1) and $peptide_length >= abs($current_var_pos2)){      
                    my @vardecoy_peptides = split(//,$vardecoy_peptide);
                    ($vardecoy_peptides[$current_var_pos2],$vardecoy_peptides[$i]) = ($vardecoy_peptides[$i],$vardecoy_peptides[$current_var_pos2]);
                    $vardecoy_peptide = join('',@vardecoy_peptides);
                    my @modrev = split /:/,$vardecoy_modif_string,-1;
                    ($modrev[$i+1],$modrev[$current_var_pos2+1]) = ($modrev[$current_var_pos2+1],$modrev[$i+1]);                      
                    $vardecoy_modif_string = join ':',@modrev;
                }
                else{
                    next;
                } 
                
                if ($peptide eq $vardecoy_peptide) {
                    next;
                }

                #calculate decoy scores for second variant
                my ($dec_masswiz_score,$dec_fractional_intensity_covered,$dec_rmsd, $dec_bycontinuity_count, $norm_var_score, $avg_var_score, $norm_var_intensity, $avg_var_intensity) = MainScoringCaller(\@exp, $vardecoy_peptide, $sum_intensity, $frag_tol, $modif_string, $charge_state, $unit, $ins_type, %variant_ions);

                push(@decoy_scores,$dec_masswiz_score);
            }
        }

        #find 2nd rank of decoy scores after normalized score
        my @sorted_decoy_scores = sort {$b <=> $a} @decoy_scores;
        if(defined $sorted_decoy_scores[0]){
            $second_rank_decoy_score = $sorted_decoy_scores[0];
            $normSecond_rank_decoy_score = $second_rank_decoy_score/((2*$peptide_length-3) + $bycontinuity_count);
        }

        #calculate the wild type score using scoring module
        ($wt_masswiz_score,$wt_fractional_intensity_covered,$wt_rmsd, $wt_bycontinuity_count,$wt_norm_var_score,$wt_avg_var_score,$wt_norm_var_intensity,$wt_avg_var_intensity) = MainScoringCaller(\@exp, $wt_peptide, $sum_intensity, $frag_tol, $modif_string, $charge_state, $unit, $ins_type, %variant_ions);

        #calculate the normalized wild type score
        $wt_normMassWizScore = $wt_masswiz_score/((2*$peptide_length-3)+$wt_bycontinuity_count);

        #normalized MassWiz Score
        my $normMassWizScore = 0;
        $normMassWizScore = $masswiz_score/((2*$peptide_length-3)+$bycontinuity_count);

        #calculate delta score
        my $deltaR2 = 0;
        my $deltaR2ratio = 0;
        if($normMassWizScore > 0 and defined $normSecond_rank_decoy_score){
            $deltaR2 = $normMassWizScore-$normSecond_rank_decoy_score;
            $deltaR2ratio = $deltaR2*100/($normMassWizScore-$deltaR2);
        }

        my $deltaMWwtRatio = 0;
        $deltaMassWizWT = 0;
        if($normMassWizScore > 0 && defined $wt_normMassWizScore ){
            $deltaMassWizWT = $normMassWizScore-$wt_normMassWizScore;
            $deltaMWwtRatio = $deltaMassWizWT*100/($normMassWizScore+$wt_normMassWizScore);
        }

        #Calculate VAS
        my $VAS = 0;
        if($normMassWizScore > 0){
            $VAS = (($normMassWizScore + $deltaR2ratio + $deltaMWwtRatio)/3)*(log($psm_count+1)/log(10))*(log($se_counts+1)/log($total_se_count));
        }
        
        #final result stored in array
        my $result_data = join("\t",$scan,$peptide,$wt_peptide,$var_type,$var_aa,$wt_aa,$var_pos,$protein_pos,$peptide_length,$protein_id,$modif_string,$charge_state,$psm_count,$se_counts,$bycontinuity_count,$fractional_intensity_covered,$masswiz_score, $normMassWizScore,$wt_normMassWizScore,$deltaR2,$deltaMassWizWT,$VAS);
        push (@results, $result_data);
        
    }
}


#to select top ranks and remove duplicates in scans based on nMW Score
my %all_results; #hash to store top ranks and adjust VAS scale
my @final_results; #2D array to store final results

my $normMassWizScore_col = 17;
my $VAS_col = 21;

my $normVAS = 0;
my $min_VAS = 0;

#Model the negative VAS values to normal distribution (for false variants) by mirroring on right side of distribution, i.e. take every negative also as a positive value so as to 
#model the null gausssian distribution with mean at zero, and SD calculated from the negative + mirror values
# This will be used to calculate pvalues for each Variant, 
# If the variant is not following this distribution, it has higher probability of being a true variant

my @VASdata;#for negative or zero values
my $VASCount=0;#for negative or zero values
my $sumsquare=0;

#if VAS scoring for duplicates is to be kept, use PSM as key, else use scan as key
if($keep_duplicates eq 'TRUE'){
    foreach my $line (@results){
        my @data = split(/\t/,$line);
        my $psm = $data[0]."-".$data[1];
        if(!exists $all_results{$psm}){
            $all_results{$psm}[0] = $data[$normMassWizScore_col];
            $all_results{$psm}[1] = $line;
        }
        else{
            if($data[$normMassWizScore_col] ne 'NA'){
                if($data[$normMassWizScore_col] > $all_results{$psm}[0]){
                    $all_results{$psm}[0] = $data[$normMassWizScore_col];
                    $all_results{$psm}[1] = $line;
                }
            }
        }
    }

}
else{
    foreach my $line (@results){
        my @data = split(/\t/,$line);
        if(!exists $all_results{$data[0]}){
            $all_results{$data[0]}[0] = $data[$normMassWizScore_col];
            $all_results{$data[0]}[1] = $line;
        }
        else{
            if($data[$normMassWizScore_col] ne 'NA'){
                if($data[$normMassWizScore_col] > $all_results{$data[0]}[0]){
                    $all_results{$data[0]}[0] = $data[$normMassWizScore_col];
                    $all_results{$data[0]}[1] = $line;
                }
            }
        }
    }
}


#calculate SD for VAS
foreach my $scan (keys %all_results){
    my @data = split(/\t/,$all_results{$scan}[1]);
    if ($data[$VAS_col] <= 0){
        $VASCount++;#double this later for mirror part
        $sumsquare+=2*($data[$VAS_col]*$data[$VAS_col]);#multiply by 2 for the mirror part
    }
}
my $mean = 0;
my $SD = sqrt($sumsquare/((2*$VASCount)-1));#2* for mirror part count

#print "SD: $SD \tCount:$VASCount\tSQ:$sumsquare\n";

#calculate zscore and pvalue for each variant
my $zed = Statistics::Zed->new(
   ccorr    => 0, #default
   tails    => 1,
   #precision_s => 3,
   precision_p => 5,
);
foreach my $scan (keys %all_results){
    my @data = split(/\t/,$all_results{$scan}[1]);
    my $vas = $data[$VAS_col];
    my $z_value = $data[$VAS_col]/$SD;
    my $p_value = $zed->z2p(value => $z_value, tails => 1);
    if($vas <= 0){
        $p_value = 1-$p_value;
    }
    $all_results{$scan}[2] = $z_value;
    $all_results{$scan}[3] = $p_value;

    #print "V:$vas\tZ:$z_value\tP:$p_value\n";<STDIN>;
}

my ($good,$ambiguous,$bad) = (0,0,0);
foreach my $key (sort keys %all_results){
    my $line = $all_results{$key}[1];
    my @data = split(/\t/,$line);

    my $VAS = $data[$VAS_col];
    my $VAS_zscore = $all_results{$key}[2];
    my $VAS_pvalue = $all_results{$key}[3];

    if($VAS_pvalue > 0.1)
    {
        #pvalue greater than 0.1 means VAS score following normal distribution 
        $data[$VAS_col+1] = $VAS_zscore;
        $data[$VAS_col+2] = $VAS_pvalue;
        $data[$VAS_col+3] = "Doubtful";
        $bad++;
    }
    elsif($VAS_pvalue > 0.01)
    {
        #when pvalue lie in the range of 0 to 0.1, its in ambiguous zone
        $data[$VAS_col+1] = $VAS_zscore;
        $data[$VAS_col+2] = $VAS_pvalue;
        $data[$VAS_col+3] = "Semi-Confident";
        $ambiguous++;
    }
    else
    {
        #when pvalue is zero, it means that the variant is not following normal distribution, hence it is a confident true variant
        $data[$VAS_col+1] = $VAS_zscore;
        $data[$VAS_col+2] = $VAS_pvalue;
        $data[$VAS_col+3] = "Confident";
        $good++;
    }

    #store in final results
    push(@final_results,\@data);
  
}

#----------------------------------
# Isobaric check and annotation
#----------------------------------
# calling Isobaric check
my $final_results_ref = Isobaric_Check(\@final_results, $frag_tol);


#printing results to file
open (OUT, ">",$OutputFile) or die $!;

print OUT "Scan\tVariantPeptide\tWildtypePeptide\tVariants\tVarAA\tWtAA\tVarPos\tProteinPos\tPeptideLength\tProteinID\tModifications\tCharge\tPSM Counts\tSE Counts\t".
"by-Continuity Count\tFractional Intensity Covered\tMassWiz Score\tNormalized MassWizScore\tNormalized WT_MassWizScore\tDeltaR2\tDeltaMassWizWT\tVAS\tZscore\tP-value\tVAS Annotation\tVariant_SubType\tVariant_Type\n";

foreach my $line (@$final_results_ref){
    print OUT join("\t",@$line),"\n";
}
close OUT;

#----------------------------------
# Disease Annotation
#----------------------------------
# calling BioAnnotation
if ($BioAnno_check eq 'TRUE'){
    my $BioAnnoOutput = $OutputFile;
    $BioAnnoOutput =~ s/\.tsv$/_BioAnno.tsv/;
    BioAnno($OutputFile, $BioAnno_file, $BioAnnoOutput);
}

#Job finished
my $t1=time();
print "Variant Scoring completed in ",sprintf("%.2f",($t1-$t0)/60),"min.\n\n";
print "Total scans rescored:".(keys %all_results).". $k scans skipped for insufficient peaks.\n";
print "$good Confident, $ambiguous Semi-Confident and $bad Doubtful variants are identified.\n";
print "Results are stored in $OutputFile\n";

print LOG "$good Confident, $ambiguous Semi-Confident and $bad Doubtful variants are identified.\n";
print LOG "Variant Scoring completed in ",sprintf("%.2f",($t1-$t0)/60),"min.\n\n";
print LOG "Results are stored in $OutputFile\n\n";
print "\t","="x60,"\n\n";

exit;