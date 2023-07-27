#!/usr/bin/perl
use strict;
use warnings;

my $usgae = "Usage: perl Plabel_batch_spectra.pl <MassWiz file> <MGF file> <output file>\n";

my $input=$ARGV[0] or die $usgae; #PgxSAVy output file
my $out=$input;
$out=~s/\.tsv/\.plabel/;

my $mgf=$ARGV[1]; #Full path to mgf file

my @mods=( 'Cys_CAM' , 'Oxidation_M'); #modifications used in the search,


#define columns in the input file
my $scan_col = 0;
my $peptide_col = 1;
my $modifications_col = 10;
my $score_col = 21;

#generally used mods in plabel format, more can be added later, as per reqmnt, This is a global list
my %modifications=(
											'Cys_CAM'=>'Carbamidomethyl[C]', #MWcode => Plabel code
											#'DEAM_N'=>'_Deamidation_NQ_0.984',
											#'DEAM_Q'=>'_Deamidation_NQ_0.984',
											#'PYRR'=>'_Gln->pyro-Glu_-17.026',
											'Oxidation_M'=>'Oxidation[M]',
											#'ACET_nterm'=>'Acetyl[AnyN-term]'
										);
 

my $header='[FilePath]'."\n"."File_Path=".$mgf."\n".'[Modification]'."\n";
my $i=1;
foreach (@mods)
	{
		$header.=$i.'='.$modifications{$_}."\n";
		$i++;
	}
$header.='[xlink]'."\n".'xlink=NULL'."\n".'[Total]'."\n".'total=100'."\n"; #TODO  add PSM count here


#############################################################################
print "Analysing input file: $input\n";

open MW,$input or die "Cannot open $input: $!";
open OUT,">$out" or die "Cannot open $out: $!";
#print $header;
my $count=0;
my $h=<MW>; #header flag

my $contents;
while(my $line = <MW>)
	{
		chomp $line;
		my @arr=split"\t",$line;
		$count++;
		#print "$line\n";
		my $spectrumID=uc($arr[$scan_col]);
		$contents.='[Spectrum'.$count.']'."\n";
		$contents.='name='.$spectrumID."\n";
		$contents.='pep1=0 '.$arr[$peptide_col].' '.$arr[$score_col];#0 PEPSEQ Score Modpos,ModCode
		my $modifs=$arr[$modifications_col];
		#print "\$mods=$mods-----";
		my $modpos;
		if ($modifs=~/[^:]/)
		{
			#print"$modifs---";
			my $temp=$modifs;
			my $c=1; #index no. of modificaton starting from 1 as defined in @mods array
			foreach my $modif (@mods)
			{
				#$modifs=~s/$modif/$modifications{$modif}/g;
				$temp=~s/$modif/$c/g;
				$c++;
			}
			my @arr1=split'',$temp;
			my $position=0;
			foreach my$ar(@arr1)
			{
				if ($ar=~m/\d/)
				{
					$modpos.=' '.$position.','.$ar;
				}
				else
				{
					$position++;
				}
			}
			$contents.=$modpos."\n";
		}
		else
		{
			$contents.="\n";
		}
	}
$header=~s/total\=100/total\=$count/m;
print OUT $header;
print OUT $contents;

print "Plabel generated: $out\n";

exit;