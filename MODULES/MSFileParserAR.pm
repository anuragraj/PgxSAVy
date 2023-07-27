#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

use constant H=>1.007825;

########################################################################
#Author: Anurag Raj
#anurag.igib@gmail.com
# VERSION 1.0.0   Date 23 Sep 2022
############################################################################
# PROGRAM TO READ MS SPECTRA FILES IN A HASH FOR ACCESSING DURING SEARCH
#########################################################################

sub spectra_parser{
	my ($RelInt,$MinPeaks,$spectra_file,$scanref)=@_;
	#print $spectra_file,"\n";
	$spectra_file=~m/\.(...)$/; #get extension
	my $ext=lc($1);
	#my @ExpSpectra;#2D array structure [0]=ScanID/SpectrumID; [1]=precursor_mz; [2]=intensity; [3]=charge [4]=RT; [5]=spectrum hashref; [6]=Mr experimental;
	my %ExpSpectra; #hash of array structure [0]=ScanID/SpectrumID; [1]=precursor_mz; [2]=intensity; [3]=charge [4]=RT; [5]=spectrum hashref; [6]=Mr experimental;
    
    if ($ext eq "mgf") 
		{
			%ExpSpectra=parse_mgf($RelInt,$MinPeaks,$spectra_file,$scanref);
		}
	else
		{
			die "Unsupported spectra file format\nOnly mgf are allowed\n";
		}
	return \%ExpSpectra;
}

###########################################################################################
###########################################################################################
sub parse_mgf{
	my ($RelInt,$MinPeaks,$spectra_file,$scanref)=@_;
	my %scans=%$scanref;
	my $flag=0;
	my $j=0; # spectrum 2-D for each scan
    my %ExpSpectra; # hash for storing spectra
	my @spectrum;
	my $flag_last=0;

	my $title="";
	open MGF,$spectra_file or die "$spectra_file not found\n";
	#open (SELECTED,">".$spectra_file."_SELECTED") or die("cannot create output file");
	while(my $line = <MGF>)
	{
		chomp $line;
		$line =~ s/\r//g;
		next if ($line =~ m/^\s+$/);
		
		if ($line =~ m/^BEGIN\sIONS$/)
		{
			$flag=1;
			next;
		}
		elsif($line =~ m/^END\sIONS$/ and $flag==1)
		{
			$flag=0;
			my @temp=@spectrum;
			#print "Before Spectrum: @temp\n";
			#if charge is missing in MGF file
			# if(!$ExpSpectra{$title}[6]){ #CHECK AGAIN
			# 	$ExpSpectra{$title}[6] = 2;
			# 	print "$title\t$ExpSpectra{$title}[6]\n";
			# }
			@temp=intensity_filter($RelInt,$MinPeaks,\@temp,$ExpSpectra{$title}[6]);
			#print "After Spectrum: @temp\n";<>;
			if ((scalar@temp) >= $MinPeaks) 
			{
				my $spec_string='';
				foreach my $value (@temp)
				{
					$spec_string.="$value->[0] $value->[1]\n";
				}
				$ExpSpectra{$title}[5]=$spec_string;
				$ExpSpectra{$title}[7]=scalar@temp; #number of peaks
				@temp=();
				$spec_string='';
				@spectrum=();
				$flag_last=0;
				$j=0;
				next;
			}
			else
			{
				$j=0;
				@spectrum=();
				$flag_last=1;
				next;
			}

		}
		######################
		#print "FLAG:$flag\nLine: '$line'\n";<STDIN>;
		if ($flag==1)	#add spectral info
		{
			if ($line =~ m/^\d+/) 
			{
				my @sp = split(/\s+/,$line);
				$spectrum[$j][0]=$sp[0];
				$spectrum[$j][1]=$sp[1];
				$j++;
			}
			elsif ($line =~ m/^TITLE/) 
			{
				#my @arr = split("=",$line);
				my $arr = $line =~ s/TITLE=//r;
				$title=$arr;
				#print "$title\n";
				if(exists $scans{$title})
				{
					$ExpSpectra{$title}[0]=();
				}
				else
				{
					$flag=0;
					#print "Scan $title not found in the input file\n";
				}
				
			}
			elsif ($line =~ m/^PEPMASS/) 
			{
				my @arr = split(/=/,$line);
				if ($arr[1]=~m/\s/g) 
				{
					my @arr2=split/\s/,$arr[1];
					$ExpSpectra{$title}[1]=$arr2[0]; #mz
					$ExpSpectra{$title}[2]=$arr2[1]; #intensity
				}
				else
				{
					$ExpSpectra{$title}[1]=$arr[1]; #mz
					$ExpSpectra{$title}[2]=undef; #intensity
				}
			}
			elsif ($line =~ m/^CHARGE/) 
			{
				my @arr = split(/=/,$line);
				$ExpSpectra{$title}[3]=$arr[1];
				$ExpSpectra{$title}[3]=~s/\+//g;
				$ExpSpectra{$title}[3]=~s/\-//g;
				$ExpSpectra{$title}[6]=($ExpSpectra{$title}[1]*$ExpSpectra{$title}[3])-($ExpSpectra{$title}[3]*H); #Exp Mr (Deconvolution step)
			}
			elsif ($line =~ m/^RTINSECONDS/) 
			{
				my @arr=split(/=/,$line);
				$ExpSpectra{$title}[4]=$arr[1];
			}
		}
		elsif($flag==0)
		{
			next;
		}

	}
	if ($flag_last==1) {
		delete($ExpSpectra{$title});
	}
	#print "$title\n";<STDIN>;
	close MGF;
	return %ExpSpectra;
}
###########################################################################################
###########################################################################################

sub intensity_filter 
	{
		my ($RelInt,$minpeaks,$ref,$precursor_mass)=@_;
		my %sp;
		foreach (@{$ref}) 
			{
				$sp{$_->[0]}=$_->[1];
			}
		my@intensity=sort{$a<=>$b}values%sp;
		my ($min,$max) = (sort {$a<=>$b}@intensity)[0,-1];
		#print"\$precursor_mass=$precursor_mass\n";
		foreach my $key (keys(%sp))
				{
					if($key>$precursor_mass-5)
					{
						delete $sp{$key};
					}
				}

	my $total_ions_required=int($precursor_mass/110)*4;
	#print "\$precursor_mass=$precursor_mass\t\$total_ions_required 1 =$total_ions_required\n";
	if(($total_ions_required%5)>0)
		{
			$total_ions_required=(int($total_ions_required/5)+1)*5;
		}
	my $bin_number=$total_ions_required/5;
	#print "\$bin_number=$bin_number\n";
	#print "\$total_ions_required=$total_ions_required\t\$bin_number=$bin_number\n\n";
	my ($min_mass,$max_mass)=(sort{$a<=>$b} keys(%sp))[0,-1];
	my $increment=($max_mass-$min_mass)/$bin_number;
	my %sp_mod;
	for(my $bin_value=$min_mass;$bin_value<$max_mass;$bin_value+=$increment)
	{
		my @sub_array=();
		my $i=0;
		foreach my $key(sort{$a<=>$b} keys(%sp)) 
			{
				if(($key>=$bin_value)&&($key<$bin_value+$increment))
				{
					$sub_array[$i][0]=$key;
					
					$sub_array[$i][1]=$sp{$key};
					$i++;
				}
				else
				{
					next;
				}
			}
		my @sorted_sub_array=sort{$b->[1]<=>$a->[1]}@sub_array;
		my $max_limit=scalar(@sorted_sub_array);
		if($max_limit>5)
			{
				$max_limit=5;
			}
		for(my $j=0;$j<$max_limit;$j++)
			{
				$sp_mod{$sorted_sub_array[$j][0]}=$sorted_sub_array[$j][1];
			}
	}
	foreach my $in (keys%sp_mod) 
		{
			if (($sp_mod{$in}*100/$max)>= $RelInt)
				{
					next;
				}
			else{
					delete $sp_mod{$in};
				}
		}		
		
	if (scalar(keys %sp_mod)>=$minpeaks)
		{
			my @spec;
			my $i=0;
			foreach (sort{$a<=>$b}keys%sp_mod) {
				$spec[$i][0]= $_;
				$spec[$i][1]= $sp_mod{$_};
				$i++;
			}
			return @spec;
		}
	else
		{
			my @spec;
			return @spec;
		}
}
1;