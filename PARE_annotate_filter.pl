#!/usr/bin/perl -w
use strict;

#Output columns of All_Reps_5cpm_Filtered_Annotated_coordinates.txt
#"Chr\tSource\tBiotype\tStart\tEnd\t\tOri\t\tPeak_Pos\tCol0_cpm1\txrn4_cpm1\tdne1xrn4_cpm1\tLog2FC_xrn4_v_Col0_Rep1\tLog2FC_xrn4_v_dne1xrn4_Rep1\t";
#"Col0_cpm2\txrn4_cpm2\tdne1xrn4_cpm2\tLog2FC_xrn4_v_Col0_Rep2\tLog2FC_xrn4_v_dne1xrn4_Rep2\tmRNA\tLocus\tAnnotation\n";
#

my $tabfile = "All_Seedling_Reps_1cpm_Filtered_Annotated.txt"; #Output of bedtools intersect that includes genomic coordinates
my $output = "All_Seedling_Reps_1cpm_Filtered_Annotated_coordinates.txt"; #Final output file with All PARE scores with new coordinats to extract Sequence


open(TABFILE, $tabfile) || die "Can't open $tabfile";
open (OUT, ">$output") || die "Can't open $output";

#Read headers
<TABFILE>;

my %hTargets = ();
my @arry = ();
my %hLocus = ();

my $loc;
my ($c1, $x1, $d1, $c2, $x2, $d2); #cpm values
my $strand;
my $ori;

############Read file into Hash tables ################

while (<TABFILE>) {
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@arry=split(/\t/, $_);
	
	if($arry[23] eq "mRNA" || $arry[23] eq "ncRNA" || $arry[23] eq "ncRNA_gene"){
		$hTargets{$arry[4]}{'chr'} = $arry[0];						
		$hTargets{$arry[4]}{'start'} = $arry[1];
		$hTargets{$arry[4]}{'end'} = $arry[2];
		$hTargets{$arry[4]}{'ori'} = $arry[3];

		$c1 = $arry[5] + $arry[6]; #Sum F + R to get a total CPM value Rep 1
		$x1 = $arry[7] + $arry[8];
		$d1 = $arry[9] + $arry[10];
						
		$c2 = $arry[13] + $arry[14]; #Sum F + R to get a total CPM value Rep 2
		$x2 = $arry[15] + $arry[16];
		$d2 = $arry[17] + $arry[18];
							
							
		$hTargets{$arry[4]}{'Col_cpm1'} = $c1;
		$hTargets{$arry[4]}{'Col_cpm2'} = $c2;
		$hTargets{$arry[4]}{'xrn4_cpm1'} = $x1;
		$hTargets{$arry[4]}{'xrn4_cpm2'} = $x2;
		$hTargets{$arry[4]}{'dbl_cpm1'} = $d1;
		$hTargets{$arry[4]}{'dbl_cpm2'} = $d2;
						
		$hTargets{$arry[4]}{'log2_CX1'} = $arry[11]; #Log2 FC xrn4 Vs Col0 Rep1 and Rep2
		$hTargets{$arry[4]}{'log2_CX2'} = $arry[19];
						
		$hTargets{$arry[4]}{'log2_MX1'} = $arry[12]; #Log2 FC xrn4 Vs dbl Rep1 and Rep2
		$hTargets{$arry[4]}{'log2_MX2'} = $arry[20];
		
		
		$hTargets{$arry[4]}{'Col_cpm1'} = $c1;
		$hTargets{$arry[4]}{'Col_cpm2'} = $c2;
		$hTargets{$arry[4]}{'xrn4_cpm1'} = $x1;
		
		$hTargets{$arry[4]}{'rChrome'} = $arry[21];
		$hTargets{$arry[4]}{'rType'} = $arry[23];
		$hTargets{$arry[4]}{'rStart'} = $arry[24];
		$hTargets{$arry[4]}{'rEnd'} = $arry[25];
		$hTargets{$arry[4]}{'rStrand'} = $arry[27];
		$hTargets{$arry[4]}{'rScore'} = $arry[28];
		$hTargets{$arry[4]}{'rDesc'} = $arry[29];
		
			
	} #End If
} #End while


foreach my $key (sort keys %hTargets){
	
	my $annot = $hTargets{$key}{'rDesc'};
	my @mRNA1 = split(/;/, $annot);
	$mRNA1[0] =~ m/^(AT\wG\d+\.\d+)/;
	my @mRNA2 = split(/:/, $mRNA1[0]);
	my $id = $mRNA2[1];
	chomp($id);
	
    #Adjust 5'P position to get 20 nt up and downstream of peak site
	my $newstart;
	my $newend;
	$loc = "$hTargets{$key}{'rChrome'}" . ":" . "$hTargets{$key}{'rStart'}" . "-" . "$hTargets{$key}{'rEnd'}";
	chomp($loc);
	#$newend = $hTargets{$key}{'start'} + 20;
	#$newstart = $hTargets{$key}{'start'} - 20;
	
	if($hTargets{$key}{'rStrand'} eq "+"){
		$newend = $hTargets{$key}{'start'} + 20;
		$newstart = $hTargets{$key}{'start'} - 19;
	} elsif ($hTargets{$key}{'rStrand'} eq "-"){
		$newend = $hTargets{$key}{'start'} + 21; #Gives 20 nt upstream sequence
		$newstart = $hTargets{$key}{'start'} - 18; #Gives 20 nt downstream sequence
	}
	
    print OUT "$hTargets{$key}{'chr'}\t" . "araport11" . "\t$hTargets{$key}{'rType'}\t$newstart\t$newend\t" . "." . "\t$hTargets{$key}{'rStrand'}\t$hTargets{$key}{'rScore'}\t$key\t";
	print OUT "$hTargets{$key}{'Col_cpm1'}\t$hTargets{$key}{'xrn4_cpm1'}\t$hTargets{$key}{'dbl_cpm1'}\t$hTargets{$key}{'log2_CX1'}\t$hTargets{$key}{'log2_MX1'}\t";
    print OUT "$hTargets{$key}{'Col_cpm2'}\t$hTargets{$key}{'xrn4_cpm2'}\t$hTargets{$key}{'dbl_cpm2'}\t$hTargets{$key}{'log2_CX2'}\t$hTargets{$key}{'log2_MX2'}\t";
	print OUT "$id\t$loc\t$annot\n";

    
}
close(OUT);

###########################################################