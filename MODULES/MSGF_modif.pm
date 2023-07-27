#!/usr/bin/perl

use strict;
use warnings;

########################################################################
#Author: Anurag Raj
#anurag.igib@gmail.com
#VERSION 1.0.0   Date 13 june 2023
############################################################################
# MODULE TO GENERATE MSGF+ MODIF STRING 
#########################################################################


sub MSGF_modif
{
	my ($peptide,$modif,$fixed_mod_str,$var_mod_str)=@_;

	my %fixed_mod_defined=();
	my %variable_mod_defined=();
	my %mod_list=ModList();
	my %pos_mod=();
	my $nterm;
	my $modif_string;
	$peptide=uc($peptide);
	
	for(my $i=0;$i<length($peptide);$i++)
		{
			$pos_mod{$i}='';
		}
	
	if($fixed_mod_str ne '')
		{
			my @fixed_mods=();	
			if(defined($fixed_mod_str))
			{		
				@fixed_mods=split(/,/,$fixed_mod_str);	
			}
			foreach my $fixed_mod_id(@fixed_mods)
				{				
					$fixed_mod_defined{$fixed_mod_id}[0]=$mod_list{$fixed_mod_id}[0];
					$fixed_mod_defined{$fixed_mod_id}[1]=$mod_list{$fixed_mod_id}[1];
					$fixed_mod_defined{$fixed_mod_id}[2]=$mod_list{$fixed_mod_id}[2];
					$fixed_mod_defined{$fixed_mod_id}[3]=$mod_list{$fixed_mod_id}[3];
					$fixed_mod_defined{$fixed_mod_id}[4]=$mod_list{$fixed_mod_id}[4];
					#print "$fixed_mod_defined{$fixed_mod_id}[0]\t$fixed_mod_defined{$fixed_mod_id}[1]\t$fixed_mod_defined{$fixed_mod_id}[2]\t$fixed_mod_defined{$fixed_mod_id}[3]\t$fixed_mod_defined{$fixed_mod_id}[4]\t$peptide","\n";
				}		
			foreach my $f_mod(values(%fixed_mod_defined))
				{						
					my $pos=-1;
					while(($pos=index($peptide,$f_mod->[4],$pos+1))>-1)
					{						
						$pos_mod{$pos}=$f_mod->[2];	
					}
				}
		}
	if($modif ne "" and $modif ne "None")
		{
			my @var_mods=();
			if(defined($var_mod_str))	
				{
					@var_mods=split(/,/,$var_mod_str);
				}
			my @amino_acids=split(//,$peptide);
			foreach my $var_mod_id(@var_mods)
				{
					$variable_mod_defined{$var_mod_id}[0]=$mod_list{$var_mod_id}[0];
					$variable_mod_defined{$var_mod_id}[1]=$mod_list{$var_mod_id}[1];
					$variable_mod_defined{$var_mod_id}[2]=$mod_list{$var_mod_id}[2];
					$variable_mod_defined{$var_mod_id}[3]=$mod_list{$var_mod_id}[3];
					$variable_mod_defined{$var_mod_id}[4]=$mod_list{$var_mod_id}[4];
					$variable_mod_defined{$var_mod_id}[5]=$mod_list{$var_mod_id}[5];
                    #print "$variable_mod_defined{$var_mod_id}[0]\t$variable_mod_defined{$var_mod_id}[1]\t$variable_mod_defined{$var_mod_id}[2]\t$variable_mod_defined{$var_mod_id}[3]\t$variable_mod_defined{$var_mod_id}[4]\t$variable_mod_defined{$var_mod_id}[5]\t$peptide","\n";
				}
			
			my @var_modifications=split(/\|/,$modif);
			foreach my $var_mod(@var_modifications)
				{	
					my @var_mod_val=split(/:/,$var_mod);
					my $pos = 0;
					if(defined $var_mod_val[1])
					{
						$pos=$var_mod_val[1]-1;
					}else{
						$nterm = $var_mod_val[0];
						next;
					}
					my $mod_desc=$var_mod_val[0];		
					foreach my $mod_ref(values(%variable_mod_defined))
					{
						if($mod_ref->[3] eq $mod_desc)
							{
								$pos_mod{$pos}=$mod_ref->[2];
							}
					}
				}	
		}
	if(defined $nterm){
		$modif_string="$nterm:";
	}else{
		$modif_string=':';
	}
	
	for(my $j=0;$j<length($peptide);$j++)
		{
			$modif_string.=$pos_mod{$j}.':';	
		}
	return($modif_string);
}


sub ModList
{
	my $mods_file='./MODULES/mods_list.txt';
	my %mods_list=();
	open(MOD,$mods_file)or die("mods_list.txt file missing in the MODULES directory");
	while(my $mod_line=<MOD>)
	{
		chomp $mod_line;
		if($mod_line=~/^\#/)
		{
			next;
		}
		else
		{
			my @line_data=split(/\t+/,$mod_line);
			my $mod_id=$line_data[1];					#MOD_ID
			$mod_id=~s/\r|\s|\n//g;
			
			my $mod_exp=$line_data[10].'@'.$line_data[4];
			$mod_exp=~s/\r|\s|\n//g;
			
			
			$mods_list{$mod_id}[0]=$line_data[2];		#Description
			$mods_list{$mod_id}[1]=$mod_exp;			#exp for TANDEM
			$mods_list{$mod_id}[2]=$line_data[5];		#MOD
			$mods_list{$mod_id}[3]=$line_data[10];		#Monoisotopic mass
			$mods_list{$mod_id}[4]=$line_data[4];		#Site
			$mods_list{$mod_id}[5]=$line_data[7];		#PSI-MS Name
		}
	}
	close MOD;
	return(%mods_list);
}

1;