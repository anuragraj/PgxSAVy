#!/usr/bin/perl

use strict;
use warnings;

########################################################################
#Author: Anurag Raj
#anurag.igib@gmail.com
#VERSION 1.0.0   Date 10 june 2023
############################################################################
# MODULE TO FILTER ISOBARIC VARIANTS AND CLASSIFICATION
#########################################################################
use constant H2O => 18.010564683;

#subroutine to check isobaric and annotate
sub Isobaric_Check{

    #Input: results array, amino acid hash, isobaric hash
    my ($final_results_ref, $mass_tolerance) = @_;

    my @final_results = @$final_results_ref;
    
    # create a modification hash
    my $modification_ref = read_modifications();
    my %modification_list = %$modification_ref;

    #amino acid table
    my %amino_acid_list=(
        #aa	Monoisotopic Mass
        "A"	=>	71.03711,
        "R"	=>	156.10111,
        "N"	=>	114.04293,
        "D"	=>	115.02694,
        "C"	=>	103.00919,
        "E"	=>	129.04259,
        "Q"	=>	128.05858,
        "G"	=>	57.02146,
        "H"	=>	137.05891,
        "I"	=>	113.08406,
        "L"	=>	113.08406,
        "K"	=>	128.09496,
        "M"	=>	131.04049,
        "F"	=>	147.06841,
        "P"	=>	97.05276,
        "S"	=>	87.03203,
        "T"	=>	101.04768,
        "W"	=>	186.07931,
        "Y"	=>	163.06333,
        "V"	=>	99.06841
    );

    #modification mass table
    my %mods = (
        "C" => {"Cys_CAM" => 57.021464},
        "M" => {"Oxidation_M" => 15.994915},
        "N" => {"Deamidation_N" => 0.984016},
        "Q" => {"Deamidation_Q" => 0.984016},
        "S" => {"PHOS_ST" => 79.966331},
        "T" => {"PHOS_ST" => 79.966331},
        "Y" => {"PHOS_STY" => 79.966331},

    );

    foreach my $psm_result (@final_results){

        my $peptide = $psm_result->[1];
        my $wt_peptide = $psm_result->[2];

        my $modif_string = $psm_result->[9];
        
        my ($wt1,$mut1,$wt2,$mut2);
        
        if($psm_result->[3] eq "Single"){
            $mut1 = $psm_result->[4];
            $wt1 = $psm_result->[5];
        }
        else{
            my @var_aa = split(/\;/,$psm_result->[4]);
            $mut1 = $var_aa[0];
            $mut2 = $var_aa[1];

            my @wt_aa = split(/\;/,$psm_result->[5]);
            $wt1 = $wt_aa[0];
            $wt2 = $wt_aa[1];
        }

        #output return
        my ($Variant_SubType, $Variant_Type) = ("NA", "NA");

        #parse modif string (::::::Oxidation_M:::Cys_CAM::Oxidation_M:::::::::)
        my @modif_array = split(/\:/,$modif_string);

        ## Creating hashes of amino acid masses and aa difference for Isobaric Filters
        my ($aa_ref, $iso_ref) = read_aa_iso_diff($peptide, \@modif_array, \%amino_acid_list, \%mods);
        my %isobaric_list = %$iso_ref;

        #calculate mass of peptide
        my $peptide_mass=peptide_mass($peptide, %amino_acid_list);
        my $wt_peptide_mass=peptide_mass($wt_peptide, %amino_acid_list);

        #check the mass difference between variant peptide and wild type peptide
        my $delta_mass= ($peptide_mass - $wt_peptide_mass);

        #print "VP:$peptide_mass\tWT:$wt_peptide_mass\tDel:$delta_mass\n";

        #check first if mutation can be justified with modification
        #Single variant
        if($psm_result->[3] eq "Single"){

            #check if delta mass or aa mass difference is within mass tolerance
            if (abs($delta_mass) <= $mass_tolerance) #Isobaric
            {
                #ISOBARIC
                #print OUTPUT "$line\t$delta_mass\tISO\tIsobaricAA\t$isobaric_list{$mut1}{$wt1}\n";
                ($Variant_SubType, $Variant_Type) = ("ISO", "IsobaricAA");
            }
            else
            {
                #check if mutation can be justified with modification
                my $matchflag=0; #if matches mod mass, it will be 1
                foreach my $modmass (keys %{ $modification_list{$wt1} })
                {
                    #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides

                    if (abs($delta_mass - $modmass) <= $mass_tolerance) 
                    {
                        #modification
                        $matchflag=1;
                        #print OUTPUT "$line\t$delta_mass\tMOD\tIsobaricMod\t$modmass\n";
                        ($Variant_SubType, $Variant_Type) = ("MOD", "IsobaricMod");
                        last;
                    }
                }
                if($matchflag==0) # no mod matched, means it is a true variant
                {
                    # it is variant
                    #print OUTPUT "$line\t$delta_mass\tVAR\tSingleVariant\tNA\n";
                    ($Variant_SubType, $Variant_Type) = ("VAR", "SingleVariant");
                    
                }
            }
            
        }
        else #two variants
        {
            if(abs($delta_mass) <= $mass_tolerance)
            { #ISO-ISO, VAR-VAR, MOD-MOD, VAR-MOD, MOD-VAR
                #if two mutations are present, check if both delta aa lies in tolerance
                if((abs($isobaric_list{$mut1}{$wt1}) <= $mass_tolerance)) # V1 = isobaric
                {
                    if ((abs($isobaric_list{$mut2}{$wt2}) <= $mass_tolerance)&&(abs($isobaric_list{$mut1}{$wt1}+$isobaric_list{$mut2}{$wt2})<=$mass_tolerance)) # 3rd condition to ensure both isobaric also under tol
                    {
                        #ISOBARIC
                        #print OUTPUT "$line\t$delta_mass\tISO-ISO\tIsobaricAA\t$isobaric_list{$mut1}{$wt1}\;$isobaric_list{$mut2}{$wt2}\n";
                        ($Variant_SubType, $Variant_Type) = ("ISO-ISO", "IsobaricAA");
                    }
                }
                else #V1 not isobaric
                {
                    #check if mutation can be justified with modification
                    my $matchflag=0; #if matches mod mass, it will be 1
                    foreach my $modmass1 (keys %{ $modification_list{$wt1}}) #V1 mod check
                    {
                        foreach my $modmass2 (keys %{ $modification_list{$wt2} }) #V2 mod check
                        {
                            #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides
                            if (abs($delta_mass - $modmass1 - $modmass2) <= $mass_tolerance) #both mod
                            {
                                #modification
                                $matchflag=1;
                                #print OUTPUT "$line\t$delta_mass\tMOD-MOD\tIsobaricMod\t$modmass1\;$modmass2\n";
                                ($Variant_SubType, $Variant_Type) = ("MOD-MOD", "IsobaricMod");
                                last;
                            }
                        }
                    }
                    if($matchflag==0) # no mod1+mod2 matched, means it is a mod-var/var-mod/var-var
                    {
                        foreach my $modmass1 (keys %{ $modification_list{$wt1} }) #MOD-VAR
                        {
                            if(abs($delta_mass - $modmass1-($isobaric_list{$mut2}{$wt2})) <= $mass_tolerance)
                            {
                                $matchflag=1;
                                #print OUTPUT "$line\t$delta_mass\tMOD-VAR\tIsobaricMod\t$modmass1\;$isobaric_list{$mut2}{$wt2}\n";
                                ($Variant_SubType, $Variant_Type) = ("MOD-VAR", "SingleVariant");
                                last;
                            }
                        }
                        
                        if($matchflag==0) # no mod1 matched, means it is a var-mod or var-var
                        {
                            foreach my $modmass2 (keys %{ $modification_list{$wt2} })# VAR MOD
                            {
                                if(abs($delta_mass -$modmass2-($isobaric_list{$mut1}{$wt1})) <= $mass_tolerance) 
                                {
                                    $matchflag=1;
                                    #print OUTPUT "$line\t$delta_mass\tVAR-MOD\tIsobaricMod\t$modmass2\;$isobaric_list{$mut1}{$wt1}\n";
                                    ($Variant_SubType, $Variant_Type) = ("VAR-MOD", "SingleVariant");
                                    last;
                                }
                            }
                            if($matchflag==0) # var-var
                            {
                                if (abs($delta_mass-$isobaric_list{$mut1}{$wt1}-$isobaric_list{$mut2}{$wt2})<=$mass_tolerance) {
                                    $matchflag=1;
                                    #print OUTPUT "$line\t$delta_mass\tVAR-VAR\tDoubleVariant\t$isobaric_list{$mut1}{$wt1}\;$isobaric_list{$mut2}{$wt2}\n";
                                    ($Variant_SubType, $Variant_Type) = ("VAR-VAR", "DoubleVariant");
                                }
                            }
                        }
                    }
                }
            }
            # Two variant with delta nonzero
            else #ISO-MOD, MOD-ISO,MOD-MOD, MOD-VAR, VAR-MOD, VAR-VAR, ISO-VAR, VAR-ISO
            {
                #if delta is not present in tolerance, check if it can be justified with modification
                #if first mutation delta aa lies in tolerance and second mutation can be justified with modification
                if(abs($isobaric_list{$mut1}{$wt1}) <= $mass_tolerance) #ISO-MOD, ISO-VAR
                {
                    my $matchflag=0; #if matches mod mass, it will be 1
                    if(exists $modification_list{$wt2})
                    {
                        foreach my $modmass2 (keys %{ $modification_list{$wt2} })
                        {
                            #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides
                            if (abs($delta_mass - $isobaric_list{$mut1}{$wt1}-$modmass2) <= $mass_tolerance) 
                            {
                                #modification
                                $matchflag=1;
                                #print OUTPUT "$line\t$delta_mass\tISO-MOD\tIsobaricMod\t$isobaric_list{$mut1}{$wt1}\;$modmass2\n";
                                ($Variant_SubType, $Variant_Type) = ("ISO-MOD", "IsobaricMod");
                                last;
                            }
                        }
                    }
                    #first mutation is in tolerance but second mutation is not mod
                    if($matchflag==0) # 
                    {
                        if (abs($delta_mass - $isobaric_list{$mut1}{$wt1}-$isobaric_list{$mut2}{$wt2}) <= $mass_tolerance)
                        {
                            #print OUTPUT "$line\t$delta_mass\tISO-VAR\tSingleVariant\t$isobaric_list{$mut1}{$wt1}\;$isobaric_list{$mut2}{$wt2}\n";
                            ($Variant_SubType, $Variant_Type) = ("ISO-VAR", "SingleVariant");
                        }
                    }
                }
                elsif(abs($isobaric_list{$mut2}{$wt2}) <= $mass_tolerance) #MOD-ISO, VAR-ISO
                {
                    #if first mutation is not in tolerance, check if second mutation is in tolerance and first mutation can be justified with modification
                    
                    my $matchflag=0; #if matches mod mass, it will be 1
                    if(exists $modification_list{$wt1})
                    {
                        foreach my $modmass1 (keys %{ $modification_list{$wt1} })
                        {
                            #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides
                            if (abs($delta_mass - $modmass1-$isobaric_list{$mut2}{$wt2}) <= $mass_tolerance)
                            {
                                #modification
                                $matchflag=1;
                                #print OUTPUT "$line\t$delta_mass\tMOD-ISO\tIsobaricMod\t$modmass1\;$isobaric_list{$mut2}{$wt2}\n";
                                ($Variant_SubType, $Variant_Type) = ("MOD-ISO", "IsobaricMod");
                                last;
                            }
                        }
                    }
                    if($matchflag==0)
                    {
                        if (abs($delta_mass - $isobaric_list{$mut1}{$wt1}-$isobaric_list{$mut2}{$wt2}) <= $mass_tolerance)
                        {
                            #print OUTPUT "$line\t$delta_mass\tVAR-ISO\tSingleVariant\t$isobaric_list{$mut1}{$wt1}\;$isobaric_list{$mut2}{$wt2}\n";
                            ($Variant_SubType, $Variant_Type) = ("VAR-ISO", "SingleVariant");

                        }
                    }
                }
                else
                {
                    #if both mutations are not in tolerance, check if both can be justified with modification
                    my $matchflag=0; #if matches mod mass, it will be 1
                    if(exists $modification_list{$wt1} && exists $modification_list{$wt2})
                    {
                        foreach my $modmass1 (keys %{ $modification_list{$wt1} })
                        {
                            foreach my $modmass2 (keys %{ $modification_list{$wt2} })
                            {
                                #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides
                                if (abs($delta_mass - $modmass1 - $modmass2) <= $mass_tolerance) 
                                {
                                    #modification
                                    $matchflag=1;
                                    #print OUTPUT "$line\t$delta_mass\tMOD-MOD\tIsobaricMod\t$modmass1\;$modmass2\n";
                                    ($Variant_SubType, $Variant_Type) = ("MOD-MOD", "IsobaricMod");
                                    last;
                                }    
                            }
                        }
                    }
                    if($matchflag==0) # Check single mod+var combinations , then double var
                    {
                        #if first mutation can be justified with modification and second mutation is varaint
                        my $matchflag_mod=0; #if matches mod mass, it will be 1
                        if(exists $modification_list{$wt1})
                        {
                            foreach my $modmass1 (keys %{ $modification_list{$wt1} })
                            {
                                #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides
                                if (abs($delta_mass - $modmass1-$isobaric_list{$mut2}{$wt2}) <= $mass_tolerance) 
                                {
                                    #modification
                                    $matchflag_mod=1;
                                    #print OUTPUT "$line\t$delta_mass\tMOD-VAR\tSingleVariant\t$modmass1\;$isobaric_list{$mut2}{$wt2}\n";
                                    ($Variant_SubType, $Variant_Type) = ("MOD-VAR", "SingleVariant");
                                    last;
                                }
                            }
                        }
                        if($matchflag_mod==0) # no mod matched, means it is a true variant
                        {
                            #if second mutation can be justified with modification and first mutation not in tolerance
                            my $matchflag_mod2=0; #if matches mod mass, it will be 1
                            if(exists $modification_list{$wt2})
                            {
                                foreach my $modmass2 (keys %{ $modification_list{$wt2} })
                                {
                                    #delta_mass = Varpep-WTpep => should be equal to modmass for nonvariant peptides
                                    if (abs($delta_mass - $isobaric_list{$mut1}{$wt1}- $modmass2) <= $mass_tolerance)
                                    {
                                        #modification
                                        $matchflag_mod2=1;
                                        #print OUTPUT "$line\t$delta_mass\tVAR-MOD\tSingleVariant\t$isobaric_list{$mut1}{$wt1}\;$modmass2\n";
                                        ($Variant_SubType, $Variant_Type) = ("VAR-MOD", "SingleVariant");
                                        last;
                                    }
                                }
                            }
                            if($matchflag_mod2==0) # no mod matched, means it is a true variant
                            {
                                # it is variant #CHECK
                                #print OUTPUT "$line\t$delta_mass\tVAR-VAR\tDoubleVariant\t$isobaric_list{$mut1}{$wt1}\;$isobaric_list{$mut2}{$wt2}\n";
                                ($Variant_SubType, $Variant_Type) = ("VAR-VAR", "DoubleVariant");
                            }
                        }
                    }
            
                }

            }

        }

        push (@{$psm_result}, $Variant_SubType, $Variant_Type);

    }

    return (\@final_results);

}

sub read_aa_iso_diff{

    my ($peptide, $modif_array_ref, $amino_acid_list_ref, $mods_ref) = @_;

    my %amino_acid_list = %{$amino_acid_list_ref};
    my @modif_array = @{$modif_array_ref};
    my %mods = %{$mods_ref};

    my @peptide_aa = split(//, $peptide);

    foreach my $i (0..$#modif_array)
    {
        if($modif_array[$i] ne ""){
            if($modif_array[$i] eq "Oxidation_M" && $peptide_aa[$i] eq "M")
            {
                $amino_acid_list{"M"} += $mods{"M"}{"Oxidation_M"};
            }
            elsif($modif_array[$i] eq "Cys_CAM" && $peptide_aa[$i] eq "C")
            {
                $amino_acid_list{"C"} += $mods{"C"}{"Cys_CAM"};
            }
            elsif($modif_array[$i] eq "DEAM_N" && $peptide_aa[$i] eq "N")
            {
                $amino_acid_list{"N"} += $mods{"N"}{"DEAM_N"};
            }
            elsif($modif_array[$i] eq "DEAMID" && $peptide_aa[$i] eq "N")
            {
                $amino_acid_list{"N"} += $mods{"N"}{"DEAMID"};
            }
            elsif($modif_array[$i] eq "PHOS_ST)" && $peptide_aa[$i] eq "S")
            {
                $amino_acid_list{"S"} += $mods{"S"}{"PHOS_STY"};
            }
            elsif($modif_array[$i] eq "PHOS_ST)" && $peptide_aa[$i] eq "T")
            {
                $amino_acid_list{"T"} += $mods{"T"}{"PHOS_STY"};
            }
            elsif($modif_array[$i] eq "PHOS_STY)" && $peptide_aa[$i] eq "Y")
            {
                $amino_acid_list{"Y"} += $mods{"Y"}{"PHOS_STY"};
            }
        }
    }

    #create a matrix of two amino acid mass difference from hash
    my %isobaric_list;
    foreach my $aa (keys %amino_acid_list)
    {
        foreach my $aa2 (keys %amino_acid_list)
        {
            my $mass_diff = ($amino_acid_list{$aa} - $amino_acid_list{$aa2});
            $isobaric_list{$aa}{$aa2} = $mass_diff;
        }
    }

    return (\%amino_acid_list, \%isobaric_list);
}

#read modification list
sub read_modifications{

    #read modification list and save in hash
    my %modification_list;
    open (MOD, "./MODULES/modifications.tsv") or die "Cannot open modifications.tsv in isobaric_check module\n";
    while (my $line = <MOD>)
    {
        chomp $line;
        $line =~ s/\r//g;
        my @data = split(/\t/,$line);
        my $mod_aa = $data[0];
        my $mod_mass = $data[2];
        $modification_list{$mod_aa}{$mod_mass} = ();
    }
    close MOD;

    return \%modification_list;
}



#calculate mass of peptide based on amino acid list
sub peptide_mass
{
    my ($peptide, %amino_acid_list) = @_;
    my $peptide_mass=0;
    my @peptide=split(//,$peptide);
    foreach my $aa (@peptide)
    {
        $peptide_mass+=$amino_acid_list{$aa};
    }
    $peptide_mass+=H2O;
    return $peptide_mass;    
}


=head1 Isobaric List
tag	AA1	AA2	description	mass difference (Da)
I/L	I	L	isoleucine and leucine substitution	0
N/D	N	D	deamidated asparagine is equivalent to aspartate	0
Q/E	Q	E	deamidated glutamine is equivalent to glutamate	0.000001
Q/K	Q	K	glutamine and lysine substitutionb	0.03638
M/F	M	F	oxidized methionine is equivalent to phenylalanineb	0.033014
D/E	D	E	methylated aspartate is equivalent to glutamate	0
S/T	S	T	methylated serine is equivalent to threonine	0.000001
N/Q	N	Q	methylated asparagine is equivalent to glutamine	0.000001
C/M	C	M	dimethylated cysteine is equivalent to methionine	0
P/E	P	E	dioxidation proline is equivalent to glutamate	0
S/E	S	E	acetylation serine is equivalent to glutamate	0

=head modifications
AA	Modification	Monoisotopic mass
N	Deamidated	0.984016
Q	Deamidated	0.984016
P	Dioxidation	31.989829
S	Acetylation	42.010565
S	Methylation	14.01565
N	Methylation	14.01565
M	Oxidation	15.994915
D	Methylation	14.01565
C	di-Methylation	28.0313

#isobraic list
my %isobaric_list= (
    'C' => {'M' => 0},	
    'D' => {'N' => 0,
            'E' => 0},	
    'E' => {'Q' => 0.000001,
            'D' => 0,
            'S' => 0,
            'P' => 0},	
    'F' => {'M' => 0.033014},	                     
    'I' => {'L' => 0},	
    'L' => {'I' => 0},	
    'K' => {'Q' => 0.03638},	
    'M' => {'C' => 0,
            'F' => 0.033014},
    'N' => {'D' => 0,
            'Q' => 0.000001},
    'P' => {'E' => 0},
    'Q' => {'E' => 0.000001,
            'K' => 0.03638,
            'N' => 0.000001},
    'S' => {'T' => 0.000001,
            'E' => 0},	
    'T' => {'S' => 0.000001}
);

=cut

1;