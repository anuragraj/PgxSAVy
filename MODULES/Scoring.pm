#!/usr/bin/perl
#################################################
#29 april 2010
#Author: Dhirendra Kumar
#sub msms_match_FT_ICR
#sub msms_match_maldi
#sub binary search masses
#sub msms_match_default
#sub msms_match_ESI_4SECTOR
#sub msms_match_maldi_TOF_PSD

#update: Anurag Raj
#19 jan 2023

#################################################
use strict;
use warnings;

use InSilicoSpectro;
use InSilicoSpectro::InSilico::ModRes;
use InSilicoSpectro::InSilico::MassCalculator;
InSilicoSpectro::init('InSilicoSpectro/config/insilicodef.xml');

sub MainScoringCaller
{
	my($spectrum_ref,$peptide,$sum_intensity,$msms_tol,$modif,$charge_state,$msms_tol_unit,$instrument_type, %variant_ions)=@_;

	my ($normalized_score,$fractional_intensity_covered,$rmsd,$bycontinuity_count,$var_score, $var_instensity,$var_ions_matched,$var_ions_matched_check, %var_scores);
	my ($norm_var_score,$avg_var_score,$norm_var_intensity,$avg_var_intensity)=(0,0,0,0);

	my $ions_count = keys %variant_ions;
	if($ions_count > 0){
		
		($normalized_score,$fractional_intensity_covered,$rmsd,$bycontinuity_count,%var_scores) = msms_match_ESI_VAR($spectrum_ref, $peptide, $sum_intensity, $msms_tol, $modif, $charge_state, $msms_tol_unit, %variant_ions); 

		my $var_ions_matched = keys %var_scores;
		if($var_ions_matched > 0){

			foreach my $var_ion (keys %var_scores){
				#sum all the scores
				$var_score += $var_scores{$var_ion}[0];
				$var_instensity += $var_scores{$var_ion}[1];
			}

			$norm_var_score = $var_score/$var_ions_matched;
			$avg_var_score = $var_score/$ions_count;
			$norm_var_intensity = $var_instensity/$var_ions_matched;
			$avg_var_intensity = $var_instensity/$ions_count;

			#print "Var ions: $ions_count\t$var_ions_matched\t$var_score\n";
		}
		
		
	}else{
		($normalized_score,$fractional_intensity_covered,$rmsd,$bycontinuity_count) = msms_match_ESI($spectrum_ref, $peptide, $sum_intensity, $msms_tol, $modif, $charge_state, $msms_tol_unit); 	
	}
	
	return ($normalized_score,$fractional_intensity_covered,$rmsd,$bycontinuity_count,$norm_var_score,$avg_var_score,$norm_var_intensity,$avg_var_intensity);
}


######################################################################
## msms_match_ESI
######################################################################
sub msms_match_ESI
{
	my ($spectrum_ref,$peptide, $sum_intensity, $msms_tol, $modif, $charge_state, $msms_tol_unit)=@_;
	my @exp=@$spectrum_ref;

	my %spectrum;
	if (!($modif))
	{
		$modif='';
	}else{
		chomp $modif;
	}

	my @ion_series=();
	if($charge_state > 1){
		@ion_series=('b','y','b++','y++','y-H2O','b-H2O','y-NH3','b-NH3');
	}else{
		@ion_series=('b','y','y-H2O','b-H2O','y-NH3','b-NH3');
	}
	getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>[@ion_series], spectrum=>\%spectrum);

	#matching of theoretical and experimental spectrum   
	###################################################################
	if($msms_tol_unit eq 'Da')
	{
		matchSpectrumGreedy(spectrum=>\%spectrum, expSpectrum=>\@exp, order=>['y','b','y++','b++','y-H2O','b-H2O','y-NH3','b-NH3'],minTol=>$msms_tol);
	}
	elsif($msms_tol_unit eq 'ppm')
	{
		matchSpectrumGreedy(spectrum=>\%spectrum, expSpectrum=>\@exp, order=>['y','b','y++','b++','y-H2O','b-H2O','y-NH3','b-NH3'],tol=>$msms_tol);
	}
	else
	{
		die ("FRAGMENT TOLERANCE UNIT not defined properly\n");
	}

	my $len = length($peptide);
	# N-/C-terminal fragments
	my $score=0;
	my $covered_intensity=0;
	my $b_y_count=0;
	my $bycontinuity_count=0;

	##########################################################################
	#score generation for each peptide based on which are the peaks matched
	##########################################################################
	my $delta_fragment_mass=0;
	my $rmsd = 0; #root mean square deviation
	foreach my $frag (keys(%{$spectrum{mass}{term}}))
	{	
		if($frag eq 'y')
		{
			for (my $i = 0; $i < @{$spectrum{ionType}{$frag}}; $i++)
			{
				for (my $j = $i*$len; $j < ($i+1)*$len; $j++)
				{
					if (defined($spectrum{match}{term}{$frag}[$j]))
					{
						$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
						$b_y_count++;
						$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
						$score=$score+(100/exp($delta_fragment_mass));
						if (defined($spectrum{match}{term}{$frag}[$j-1]))
							{
								$score=$score+(50/exp($delta_fragment_mass));#Bonus for continuity of y series
								$bycontinuity_count++;
							}
						else 
							{
								$score=$score-(50/exp($delta_fragment_mass));#Bonus for continuity of y series
								$bycontinuity_count++;
							}
						if (defined($spectrum{match}{term}{'y-H2O'}[$j]))
							{
								$delta_fragment_mass=abs($spectrum{match}{term}{'y-H2O'}[$j][0]-$spectrum{mass}{term}{'y-H2O'}[$j]);
								$covered_intensity+=$spectrum{match}{term}{'y-H2O'}[$j][1];
								$score=$score+(25/exp($delta_fragment_mass));
							}
						if (defined($spectrum{match}{term}{'y-NH3'}[$j]))
							{
								$delta_fragment_mass=abs($spectrum{match}{term}{'y-NH3'}[$j][0]-$spectrum{mass}{term}{'y-NH3'}[$j]);
								$covered_intensity+=$spectrum{match}{term}{'y-NH3'}[$j][1];
								$score=$score+(25/exp($delta_fragment_mass));
							}
						#calculate the root mean square deviation
						$rmsd=$rmsd+($delta_fragment_mass*$delta_fragment_mass);
					}
				}
			}
		}
		if($frag eq 'b')
		{
			for (my $i = 0; $i < @{$spectrum{ionType}{$frag}}; $i++)
			{
				for (my $j = $i*$len; $j < ($i+1)*$len; $j++)
				{						
					if (defined($spectrum{match}{term}{$frag}[$j]))
					{
						$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
						#print $delta_fragment_mass."\n";
						$b_y_count++;
						$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
						$score=$score+(100/exp($delta_fragment_mass));
						if (defined($spectrum{match}{term}{$frag}[$j-1]))
						{
							$score=$score+(20/exp($delta_fragment_mass));#Bonus for continuity for b series
							$bycontinuity_count++;
						}
						else 
						{
							$score=$score-(20/exp($delta_fragment_mass));#Bonus for continuity of y series
							$bycontinuity_count++;
						}
						if (defined($spectrum{match}{term}{'b-H2O'}[$j]))
						{
							$delta_fragment_mass=abs($spectrum{match}{term}{'b-H2O'}[$j][0]-$spectrum{mass}{term}{'b-H2O'}[$j]);
							$covered_intensity+=$spectrum{match}{term}{'b-H2O'}[$j][1];
							$score=$score+(25/exp($delta_fragment_mass));
						}
						if (defined($spectrum{match}{term}{'b-NH3'}[$j]))
						{
							$delta_fragment_mass=abs($spectrum{match}{term}{'b-NH3'}[$j][0]-$spectrum{mass}{term}{'b-NH3'}[$j]);
							$covered_intensity+=$spectrum{match}{term}{'b-NH3'}[$j][1];
							$score=$score+(25/exp($delta_fragment_mass));
						}
						#calculate the root mean square deviation
						$rmsd=$rmsd+($delta_fragment_mass*$delta_fragment_mass);
					}
				}

			}
		}
		if($frag eq 'b++')
		{	
			for (my $j = 0; $j < $len; $j++)
			{
				if (defined($spectrum{match}{term}{$frag}[$j]))
				{
					$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
					$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
					$score=$score+(100/exp($delta_fragment_mass));
					#calculate the root mean square deviation
					$rmsd = $rmsd + ($delta_fragment_mass*$delta_fragment_mass);
					$b_y_count++;
				}
			}
		}
		if($frag eq 'y++')
		{	
			for (my $j = 0; $j < $len; $j++)
			{
				if (defined($spectrum{match}{term}{$frag}[$j]))
				{
					$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
					$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
					$score=$score+(100/exp($delta_fragment_mass));
					#calculate the root mean square deviation
					$rmsd = $rmsd + ($delta_fragment_mass*$delta_fragment_mass);
					$b_y_count++;
				}
			}
		}
	}	

	my $final_score=$score*sqrt($covered_intensity/$sum_intensity);
	my $fractional_intensity_covered = $covered_intensity/$sum_intensity;
	#calculate the root mean square deviation
	if($b_y_count>0)
	{
		$rmsd=sqrt($rmsd/$b_y_count);
	}
	else
	{
		$rmsd=0;
	}
	#push(@normalized_score,$final_score, $covered_intensity, $b_y_count);
	#return \@normalized_score ;
	return ($final_score,$fractional_intensity_covered,$rmsd,$bycontinuity_count);
}


######################################################################
## msms_match_ESI_VAR
######################################################################
sub msms_match_ESI_VAR
{
	my ($spectrum_ref,$peptide, $sum_intensity, $msms_tol, $modif, $charge_state, $msms_tol_unit, %variant_ions)=@_;
	my @exp=@$spectrum_ref;

	my %spectrum;
	if (!($modif))
	{
		$modif='';
	}else{
		chomp $modif;
	}

	my @ion_series=();
	if($charge_state > 1){
		@ion_series=('b','y','b++','y++','y-H2O','b-H2O','y-NH3','b-NH3');
	}else{
		@ion_series=('b','y','y-H2O','b-H2O','y-NH3','b-NH3');
	}
	
	getFragmentMasses(pept=>$peptide, modif=>$modif, fragTypes=>[@ion_series], spectrum=>\%spectrum);

	#matching of theoretical and experimental spectrum   
	###################################################################
	if($msms_tol_unit eq 'Da')
	{
		matchSpectrumGreedy(spectrum=>\%spectrum, expSpectrum=>\@exp, order=>['y','b','y++','b++','y-H2O','b-H2O','y-NH3','b-NH3'],minTol=>$msms_tol);
	}
	elsif($msms_tol_unit eq 'ppm')
	{
		matchSpectrumGreedy(spectrum=>\%spectrum, expSpectrum=>\@exp, order=>['y','b','y++','b++','y-H2O','b-H2O','y-NH3','b-NH3'],tol=>$msms_tol);
	}
	else
	{
		die ("FRAGMENT TOLERANCE UNIT not defined properly\n");
	}

	my $len = length($peptide);
	# N-/C-terminal fragments
	my $score=0;
	my $covered_intensity=0;
	my $b_y_count=0;
	my $bycontinuity_count=0;
	my $var_score = 0;
	my $var_instensity = 0;
	my $var_ions_matched = 0;
	my %var_scores;

	##########################################################################
	#score generation for each peptide based on which are the peaks matched
	##########################################################################
	my $delta_fragment_mass=0;
	my $rmsd = 0; #root mean square deviation
	foreach my $frag (keys(%{$spectrum{mass}{term}}))
	{	
		if($frag eq 'y')
		{
			for (my $i = 0; $i < @{$spectrum{ionType}{$frag}}; $i++)
			{
				for (my $j = $i*$len; $j < ($i+1)*$len; $j++)
				{
					if (defined($spectrum{match}{term}{$frag}[$j]))
					{
						$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
						$b_y_count++;
						$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
						$score=$score+(100/exp($delta_fragment_mass));
						if (defined($spectrum{match}{term}{$frag}[$j-1]))
							{
								$score=$score+(50/exp($delta_fragment_mass));#Bonus for continuity of y series
								$bycontinuity_count++;
							}
						else 
							{
								$score=$score-(50/exp($delta_fragment_mass));#Bonus for continuity of y series
								$bycontinuity_count++;
							}
						if (defined($spectrum{match}{term}{'y-H2O'}[$j]))
							{
								$delta_fragment_mass=abs($spectrum{match}{term}{'y-H2O'}[$j][0]-$spectrum{mass}{term}{'y-H2O'}[$j]);
								$covered_intensity+=$spectrum{match}{term}{'y-H2O'}[$j][1];
								$score=$score+(25/exp($delta_fragment_mass));
							}
						if (defined($spectrum{match}{term}{'y-NH3'}[$j]))
							{
								$delta_fragment_mass=abs($spectrum{match}{term}{'y-NH3'}[$j][0]-$spectrum{mass}{term}{'y-NH3'}[$j]);
								$covered_intensity+=$spectrum{match}{term}{'y-NH3'}[$j][1];
								$score=$score+(25/exp($delta_fragment_mass));
							}
						#calculate the root mean square deviation
						$rmsd=$rmsd+($delta_fragment_mass*$delta_fragment_mass);
					}
					if(exists $variant_ions{"y".$j})
					{
						if (defined($spectrum{match}{term}{$frag}[$j-1]))
						{
							$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j-1][0]-$spectrum{mass}{term}{$frag}[$j-1]);
							$var_score=100/exp($delta_fragment_mass);
							$var_instensity = $spectrum{match}{term}{$frag}[$j-1][1];
							
							if(not exists $var_scores{'y'.$j})
							{
								$var_scores{'y'.$j}[0] = $var_score;
								$var_scores{'y'.$j}[1] = $var_instensity;
							}else{
								if($var_score > $var_scores{'y'.$j}[0])
								{
									$var_scores{'y'.$j}[0] = $var_score;
									$var_scores{'y'.$j}[1] = $var_instensity;
								}
							}
							#print "Variant y ion matched: $variant_ions{'y'.$j} $spectrum{match}{term}{$frag}[$j-1][0] $spectrum{mass}{term}{$frag}[$j-1] $var_score $var_instensity\n";
						}
					}
				}
			}
		}
		if($frag eq 'b')
		{
			for (my $i = 0; $i < @{$spectrum{ionType}{$frag}}; $i++)
			{
				for (my $j = $i*$len; $j < ($i+1)*$len; $j++)
				{						
					if (defined($spectrum{match}{term}{$frag}[$j]))
					{
						$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
						#print $delta_fragment_mass."\n";
						$b_y_count++;
						$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
						$score=$score+(100/exp($delta_fragment_mass));
						if (defined($spectrum{match}{term}{$frag}[$j-1]))
						{
							$score=$score+(20/exp($delta_fragment_mass));#Bonus for continuity for b series
							$bycontinuity_count++;
						}
						else 
						{
							$score=$score-(20/exp($delta_fragment_mass));#Bonus for continuity of y series
							$bycontinuity_count++;
						}
						if (defined($spectrum{match}{term}{'b-H2O'}[$j]))
						{
							$delta_fragment_mass=abs($spectrum{match}{term}{'b-H2O'}[$j][0]-$spectrum{mass}{term}{'b-H2O'}[$j]);
							$covered_intensity+=$spectrum{match}{term}{'b-H2O'}[$j][1];
							$score=$score+(25/exp($delta_fragment_mass));
						}
						if (defined($spectrum{match}{term}{'b-NH3'}[$j]))
						{
							$delta_fragment_mass=abs($spectrum{match}{term}{'b-NH3'}[$j][0]-$spectrum{mass}{term}{'b-NH3'}[$j]);
							$covered_intensity+=$spectrum{match}{term}{'b-NH3'}[$j][1];
							$score=$score+(25/exp($delta_fragment_mass));
						}
						#calculate the root mean square deviation
						$rmsd=$rmsd+($delta_fragment_mass*$delta_fragment_mass);
					}
					if(exists $variant_ions{"b".$j})
					{
						if (defined($spectrum{match}{term}{$frag}[$j-1]))
						{
							$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j-1][0]-$spectrum{mass}{term}{$frag}[$j-1]);
							$var_score=100/exp($delta_fragment_mass);
							$var_instensity = $spectrum{match}{term}{$frag}[$j-1][1];
							if(not exists $var_scores{'b'.$j})
							{
								$var_scores{'b'.$j}[0] = $var_score;
								$var_scores{'b'.$j}[1] = $var_instensity;
							}else{
								if($var_score > $var_scores{'b'.$j}[0])
								{
									$var_scores{'b'.$j}[0] = $var_score;
									$var_scores{'b'.$j}[1] = $var_instensity;
								}
							}
							#print "Variant b ion matched: $variant_ions{'b'.$j} $spectrum{match}{term}{$frag}[$j-1][0] $spectrum{mass}{term}{$frag}[$j-1] $var_score $var_instensity\n";
						}
					}
				}

			}
		}
		if($frag eq 'b++')
		{	
			for (my $j = 0; $j < $len; $j++)
			{
				if (defined($spectrum{match}{term}{$frag}[$j]))
				{
					$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
					$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
					$score=$score+(100/exp($delta_fragment_mass));
					#calculate the root mean square deviation
					$rmsd = $rmsd + ($delta_fragment_mass*$delta_fragment_mass);
					$b_y_count++;
				}
				if(exists $variant_ions{"b".$j})
				{
					if (defined($spectrum{match}{term}{$frag}[$j-1]))
					{
						$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j-1][0]-$spectrum{mass}{term}{$frag}[$j-1]);
						$var_score=100/exp($delta_fragment_mass);
						$var_instensity = $spectrum{match}{term}{$frag}[$j-1][1];
						if(not exists $var_scores{'b'.$j})
						{
							$var_scores{'b'.$j}[0] = $var_score;
							$var_scores{'b'.$j}[1] = $var_instensity;
						}else{
							if($var_score > $var_scores{'b'.$j}[0])
							{
								$var_scores{'b'.$j}[0] = $var_score;
								$var_scores{'b'.$j}[1] = $var_instensity;
							}
						}
						#print "Variant b++ ion matched: $variant_ions{'b'.$j} $spectrum{match}{term}{$frag}[$j-1][0] $spectrum{mass}{term}{$frag}[$j-1] $var_score $var_instensity\n";
					}
				}
			}
		}
		if($frag eq 'y++')
		{	
			for (my $j = 0; $j < $len; $j++)
			{
				if (defined($spectrum{match}{term}{$frag}[$j]))
				{
					$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j][0]-$spectrum{mass}{term}{$frag}[$j]);
					$covered_intensity+=$spectrum{match}{term}{$frag}[$j][1];
					$score=$score+(100/exp($delta_fragment_mass));
					#calculate the root mean square deviation
					$rmsd = $rmsd + ($delta_fragment_mass*$delta_fragment_mass);
					$b_y_count++;
				}
				if(exists $variant_ions{"y".$j})
				{
					if (defined($spectrum{match}{term}{$frag}[$j-1]))
					{
						$delta_fragment_mass=abs($spectrum{match}{term}{$frag}[$j-1][0]-$spectrum{mass}{term}{$frag}[$j-1]);
						$var_score=100/exp($delta_fragment_mass);
						$var_instensity = $spectrum{match}{term}{$frag}[$j-1][1];
						if(not exists $var_scores{'y'.$j})
						{
							$var_scores{'y'.$j}[0] = $var_score;
							$var_scores{'y'.$j}[1] = $var_instensity;
						}else{
							if($var_score > $var_scores{'y'.$j}[0])
							{
								$var_scores{'y'.$j}[0] = $var_score;
								$var_scores{'y'.$j}[1] = $var_instensity;
							}
						}
						#print "Variant y++ ion matched: $variant_ions{'y'.$j} $spectrum{match}{term}{$frag}[$j-1][0] $spectrum{mass}{term}{$frag}[$j-1] $var_score $var_instensity\n";
					}
				}
			}
		}
	}	

	my $final_score=$score*sqrt($covered_intensity/$sum_intensity);
	my $fractional_intensity_covered = $covered_intensity/$sum_intensity;
	#calculate the root mean square deviation
	if($b_y_count>0)
	{
		$rmsd=sqrt($rmsd/$b_y_count);
	}
	else
	{
		$rmsd=0;
	}
	#push(@normalized_score,$final_score, $covered_intensity, $b_y_count);
	#return \@normalized_score ;
	return ($final_score,$fractional_intensity_covered,$rmsd,$bycontinuity_count,%var_scores);
}



######################################################################
## binary_search_masses
######################################################################
sub binary_search_masses
{
####################################################################################################
####### Subroutine to search pepptide mass range in sorted array of precursor masses . takes care of multiple precursor masses in a range
####### returns array of positions in array#########################################################################
#####################################################################################################
	my ($x, $a,$tol,$unit) = @_;            # search for x in array a
	my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
	my $i; # index of probe
	my $lower_value;
	my $higher_value;
	if ($unit eq 'Da') 
		{
				$lower_value=$x-$tol;
				$higher_value=$x+$tol;
		}
	elsif($unit eq 'ppm')
	{
			($lower_value,$higher_value)=ppm2da($x,$tol);
	}
	else{
		die"Error in defining tolerance value\n";
	}
   	my @masses=@$a;
	my $flag=0;
	my @pos=();
    while ($l <= $u) 
		{
						if($flag==0)
						{
							$i = int(($l + $u)/2);
							if ($masses[$i][6] < $lower_value)
								{
									$l = $i+1;
								}
							elsif ($masses[$i][6] > $higher_value)
								{
									$u = $i-1;
								} 
							else
								{
									push(@pos,$i); # found
									$flag=1;
								}
						}
						else
						{
							last;
						}
	   }
	if($flag==1)
		{
			my $lower_bound=0;
			my $upper_bound=0;
			my $st=$i;
					while($lower_bound==0)
					{
						$i=$i-1;
							if (($masses[$i][6] < $lower_value)||($i<0))
								{
									$lower_bound=1;
								}
							else
								{
									push(@pos,$i); # found
								}				
					}
			my $j=$st;
					while($upper_bound==0)
							{
								$j=$j+1;
									if (($j>=$#masses)||($masses[$j][6] > $higher_value))
										{
											$upper_bound=1;
										}
									else
										{
											push(@pos,$j); # found
										}				
							}
		}
	return @pos;
}


sub max_intensity
{
		my $spectrum_ref=shift;
		my @exp=@$spectrum_ref;
		my @data=sort{$a->[1]<=>$b->[1]}@exp;
		my $max_intensity=$data[-1][1];
		return $max_intensity;
}


1;
