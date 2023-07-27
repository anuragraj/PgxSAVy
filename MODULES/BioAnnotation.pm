#!/usr/bin/perl

#############################################################################################
# Program       : BioAnnotation.pm (To assign functional classifications to Variant peptides)  #                                                    #
# Author        : Anurag Raj                                                                #
# Date          : 23.06.2023                                                                #
# Version       : 1.0                                                                       #
#############################################################################################

use strict;
use warnings;

sub BioAnno{

    my ($variant_file, $BioAnno_file, $output_file) = @_;

    my $start_run = time();

    #store variant file in hash
    my %variant_hash;

    #AA code required for db file mapping
    my %aacode = (
        'A' => 'ala',
        'R' => 'arg',
        'N' => 'asn',
        'D' => 'asp',
        'B' => 'asx',
        'C' => 'cys',
        'E' => 'glu',
        'Q' => 'gln',
        'Z' => 'glx',
        'G' => 'gly',
        'H' => 'his',
        'I' => 'ile',
        'L' => 'leu',
        'K' => 'lys',
        'M' => 'met',
        'F' => 'phe',
        'P' => 'pro',
        'S' => 'ser',
        'T' => 'thr',
        'W' => 'trp',
        'Y' => 'tyr',
        'V' => 'val',
        '*' => 'ter', #stop codon
        '-' => 'del' #deletion
    );

    open (IN, "$variant_file") or die "Can't open file $variant_file\n";
    #header line
    my $header = <IN>;
    #define columns of variant input file
    my ($scan_col, $peptide_col, $protein_id_col, $variant_type_col, $varAA_col, $wtAA_col, $varpos_col, $vas_annotation_col, $isobaric_annotation_col) = (0, 1, 9, 3, 4, 5, 7, 24, 26);

    while( my $line = <IN>){
        chomp $line;
        my @line = split("\t", $line);

        my $scanID = $line[$scan_col];
        my $peptide = $line[$peptide_col];
        my $protein_id = $line[$protein_id_col];
        #run for all protein hits of scan
        #nxp:NX_Q96JE9-2(TP:33-63)(Var:T|37|A;P|40|Q)
        if($protein_id =~ /^nxp/i){
            my $prot_id = (split(/\(/, $protein_id))[0];
            $protein_id = $prot_id;
            $protein_id =~ s/nxp:NX_//;
        }

        my $variant_type = $line[$variant_type_col];
        my $varAA = $line[$varAA_col];
        my $wtAA = $line[$wtAA_col];
        my $varpos = $line[$varpos_col];
        my $vas_annotation = $line[$vas_annotation_col];
        my $isobaric_annotation = $line[$isobaric_annotation_col];

        if($variant_type =~ /Double/){
            #two variants, check each variant
            my @var_AA = split(";", $varAA);
            my @wt_AA = split(";", $wtAA);

            my @var_pos = split(";", $varpos);
            
            $variant_hash{$protein_id}{$aacode{$wt_AA[0]}}{$aacode{$var_AA[0]}}{$var_pos[0]}[0] = $scanID;
            $variant_hash{$protein_id}{$aacode{$wt_AA[0]}}{$aacode{$var_AA[0]}}{$var_pos[0]}[1] = $peptide;
            $variant_hash{$protein_id}{$aacode{$wt_AA[0]}}{$aacode{$var_AA[0]}}{$var_pos[0]}[2] = $vas_annotation;
            $variant_hash{$protein_id}{$aacode{$wt_AA[0]}}{$aacode{$var_AA[0]}}{$var_pos[0]}[3] = $isobaric_annotation;


            $variant_hash{$protein_id}{$aacode{$wt_AA[1]}}{$aacode{$var_AA[1]}}{$var_pos[1]}[0] = $scanID;
            $variant_hash{$protein_id}{$aacode{$wt_AA[1]}}{$aacode{$var_AA[1]}}{$var_pos[1]}[1] = $peptide;
            $variant_hash{$protein_id}{$aacode{$wt_AA[1]}}{$aacode{$var_AA[1]}}{$var_pos[1]}[2] = $vas_annotation;
            $variant_hash{$protein_id}{$aacode{$wt_AA[1]}}{$aacode{$var_AA[1]}}{$var_pos[1]}[3] = $isobaric_annotation;
            
        }
        else
        {   
            $variant_hash{$protein_id}{$aacode{$wtAA}}{$aacode{$varAA}}{$varpos}[0] = $scanID;
            $variant_hash{$protein_id}{$aacode{$wtAA}}{$aacode{$varAA}}{$varpos}[1] = $peptide;
            $variant_hash{$protein_id}{$aacode{$wtAA}}{$aacode{$varAA}}{$varpos}[2] = $vas_annotation;
            $variant_hash{$protein_id}{$aacode{$wtAA}}{$aacode{$varAA}}{$varpos}[3] = $isobaric_annotation;

        }
        
    }
    close IN;    
    


    #Reading database file and storing in hash
    open(DB, "$BioAnno_file") or die "Can't open $BioAnno_file. BioAnnotation not done!!\n";
    my %db;

    print "Mapping varinats for BioAnnotation...(it will take sometime)\n";

    my $db_header;
    #to check if variant is present in database
    my $flag = 0;
    my $skip = 1;

    #write results
    open(OUT, ">$output_file") or die "Can't open to write BioAnnotation output:$output_file\n";
    print OUT "ScanID\tPeptide\tProteinID\tVariant AA Change\tVAS Annottaion\tVariant Type\tGene Name\tSourceDatabaseID\tConsequence Type\tClinical Significance\tPhenotype/Disease\tPhenotype/Disease Source\tCytogenetic Band\tChromosome Coordinate\tEnsembl gene ID\tEnsembl transcript ID\tEnsembl translation ID\tEvidence\n";

    while(my $line = <DB>){

        if($line =~ /######## VARIANT INDEX ########/){
            $db_header = <DB>;
            $skip = 0;
            next;
        }
        next if $skip == 1;

        next if $line =~ /^\s/;
        next if $line =~ /^-/;
        next if $line =~ /^_/;
        next if $line =~ /^Copyrighted/;
        next if $line =~ /^Distributed/;

        chomp $line;
        $line =~ s/\r//g;

        my @line = split("\t", $line);

        #store data in hash
        my $protein_id = $line[1];
        my $var_aa_change = $line[2]; #p.Val123Ile
        
        #parse variant aa change with letters and numbers
        $var_aa_change =~ /p\.(\D+)(\d+)(\D+)/;
        my $wt_aa = lc($1);
        my $var_aa = lc($3);
        my $var_pos = $2;

        #print "$protein_id\t$wt_aa\t$var_aa\t$var_pos\n";

        #check if DB variant is in PSM variant file
        if(exists $variant_hash{$protein_id}{$wt_aa}{$var_aa}{$var_pos}){
            #print results
            print OUT "$variant_hash{$protein_id}{$wt_aa}{$var_aa}{$var_pos}[0]\t$variant_hash{$protein_id}{$wt_aa}{$var_aa}{$var_pos}[1]\t$protein_id\t$var_aa_change\t$variant_hash{$protein_id}{$wt_aa}{$var_aa}{$var_pos}[2]\t$variant_hash{$protein_id}{$wt_aa}{$var_aa}{$var_pos}[3]\t$line[0]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$line[12]\t$line[13]\n";
            $flag = 1;
        }

    }
    close(DB);

    if($flag == 0){
        print "No Annotation match found!!!\nCheck variants in results!!!\n";
    }
    close(OUT);

    my $end_run = time();
    my $run_time = $end_run - $start_run;
    print "BioAnnottaion completed in $run_time seconds\n";
}

1;