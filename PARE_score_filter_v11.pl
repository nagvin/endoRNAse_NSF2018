#!/usr/bin/perl -w
use strict;

my $tabfile = "All_PARE_Scores.bg";
my $output = "All_Seedling_1cpm_Filtered_PARE_Scores.txt";

my @arry=();

open(TABFILE, $tabfile) || die "Can't open $tabfile";
open (OUT1, ">$output") || die "Can't open $output";

#Read headers
<TABFILE>;

my $header1 = "Chr\tStart\tEnd\tOri\tPeak_Pos\tCol0_Fwd1\tCol0_Rev1\txrn4_Fwd1\txrn4_Rev1\tdne1xrn4_Fwd1\tdne1xrn4_Rev1\tLog2FC_xrn4_v_Col0_Rep1\tLog2FC_xrn4_v_dne1xrn4_Rep1\t";
my $header2 = "Col0_Fwd2\tCol0_Rev2\txrn4_Fwd2\txrn4_Rev2\tdne1xrn4_Fwd2\tdne1xrn4_Rev2\tLog2FC_xrn4_v_Col0_Rep2\tLog2FC_xrn4_v_dne1xrn4_Rep2\n";
my $destination = "$header1" . "$header2";
my $loc;
my $ori;
############Read file into Hash tables ################
while (<TABFILE>) {
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal

	@arry=split(/\t/, $_);
	if($arry[0] ne "Mt" && $arry[0] ne "Pt"){
			if($arry[5] >=1 || $arry[6] >= 1 && $arry[13] >= 1 || $arry[14]>= 1){#filter for xrn4_F and xrn4_R Reps 1 anf 2 having values > 1 CPM
				#if($arry[9] >= 1 || $arry[10] >= 1 && $arry[17] >= 1 || $arry[18] >= 1) {#filter for log2 FC xrn4/Col0 and xrn4/dne1xrn4 > 1 log2 i.e. 2-fold change in both replicates
						$loc = "$arry[0]" . ":". "$arry[1]" . "-" . "$arry[2]";
						chomp($loc);
						if($arry[5] > $arry[6] && $arry[13] > $arry[14]){
						$ori = "-";
						} elsif ($arry[5] < $arry[6] && $arry[13] < $arry[14]){
						$ori = "+";
						}
						$destination .= "$arry[0]\t$arry[1]\t$arry[2]\t$ori\t$loc\t$arry[3]\t$arry[4]\t$arry[5]\t$arry[6]\t$arry[7]\t$arry[8]\t$arry[9]\t$arry[10]\t";
						$destination .= "$arry[11]\t$arry[12]\t$arry[13]\t$arry[14]\t$arry[15]\t$arry[16]\t$arry[17]\t$arry[18]\n";
			#	} #End If
			} #End If
	} #End If
} #End while
print OUT1 $destination;
close(OUT1);