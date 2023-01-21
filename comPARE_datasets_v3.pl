#!/usr/bin/perl -w
use strict;

#Output columns of All_Reps_5cpm_Filtered_Annotated_coordinates.txt
#"Chr\tSource\tBiotype\tStart\tEnd\t\tOri\t\tPeak_Pos\tCol0_cpm1\txrn4_cpm1\tdne1xrn4_cpm1\tLog2FC_xrn4_v_Col0_Rep1\tLog2FC_xrn4_v_dne1xrn4_Rep1\t";
#"Col0_cpm2\txrn4_cpm2\tdne1xrn4_cpm2\tLog2FC_xrn4_v_Col0_Rep2\tLog2FC_xrn4_v_dne1xrn4_Rep2\tmRNA\tLocus\tAnnotation\n";
#
#
#comPARE_v3.pl:
#Include gene expression results rpkm values and fold change values Col vs xrn4, xrn4 vs. double and Col vs. double
#Includes NMD datasets
#Include Cap Datasets (C-PARE) analyzed using the ComPARE pipeline (Hurtig et al., 2021)
#
#

my $seqfile = "Seedling_Peak_Sequences.txt"; #Output of bedtools getfasta
my $tabfile = "All_Seedling_Reps_1cpm_Filtered_Annotated_coordinates.txt"; #Final output file with All PARE scores with new coordinats to extract Sequence
my $list1 = "list1.txt"; #Isoform expression: Col-0 vs. xrn4-5
my $list2 = "list2.txt"; #Isoform expression: Col-0 vs. dne1xrn4-5
my $list3 = "list3.txt"; #Isoform expression: dne1xrn4 vs. xrn4-5
my $capsites = "CPARE_Annotated_coordinates.txt"; #CPARE analysis using the New Degradome pipeline Libs ATH292 and ATH293

########## NMD Datasets #######################
my $u1 = "upf1_up.txt"; # From Degitiar et al., 2015 #All proteincoding genes
my $u3 = "upf3_up.txt"; # From Degitiar et al., 2015 #All proteincoding genes
my $ude = "UDE.txt"; #From VK Raxena Plant Cell 2020 comparing pad4 and pad4upf1 
my $sde = "SDE.txt"; #From VK Raxena Plant Cell 2020 comparing pad4 and pad4smg7
my $sde0 = "SMG7DE.txt"; #From Riha lab Gloggnitzer 2014 Cell Microbes comparing pad4 and pad4smg7
my $ude0 = "GSE41432.txt"; #GSE41432_gene_expression_RPKM.txt Wachner lab - No q-values just RPKM
my $uORF = "TAIR10_uORFs_formatted.txt"; #From TAIR10 dataset
my $peaks = "Final_Sites.txt"; #sites identified after visually looking at these Major internals and Minor internals
##############################################
my $output = "ComPARE_dne1_Targets_Leaf_Full_NMD_Finalized.txt";


open(TABFILE, $tabfile) || die "Can't open $tabfile";
open(LX, $list1) || die "Can't open $list1";
open(LD, $list2) || die "Can't open $list2";
open(LDX, $list3) || die "Can't open $list3";

open(SEQUENCES, $seqfile) || die "Can't open $seqfile";
open(CAPS, $capsites) || die "Can't open $capsites";

open(U1, $u1) || die "Can't open $u1";
open(U3, $u3) || die "Can't open $u3";
open(UDE, $ude) || die "Can't open $ude";
open(UD, $ude0) || die "Can't open $ude0";
open(SDE, $sde) || die "Can't open $sde";
open(SMGDE, $sde0) || die "Can't open $sde0";

open(SITES, $peaks) || die "Can't open $peaks";
open(uORF, $uORF) || die "Can't open $uORF";



#######################################################
open (OUT, ">$output") || die "Can't open $output";

#Read headers
#<TABFILE>;
#<SEQUENCE>;
<LX>;
<LD>;
<LDX>;


<U1>;<U1>;
<U3>;<U3>;
<UD>;<UD>;

<UDE>;<UDE>;<UDE>;<UDE>;
<SDE>;<SDE>;<SDE>;<SDE>;
<SMGDE>; <SMGDE>;

<SITES>; #read headers
<SITES>; <SITES>;
<uORF>;


my %hTargets = ();
my %hCaps = ();

my @arry = ();
my %hSeq = ();

my %hLX =();
my %hLD =();
my %hLDX =();
my %hSites = ();
my %huORF=();

my %hU1=();
my %hU3=();
my %hUDE=();
my %hSDE=();
my %hSMGDE=();
my %hUPFDE=();
my @array=();

####################### Read Final site list and TAIR10 uORF files #############
while(<SITES>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);

	$hSites{$array[1]}{'mRNA'}=$array[0];
	$hSites{$array[1]}{'peak_type'}=$array[2];
	$hSites{$array[1]}{'major_peak'}=$array[3];
	$hSites{$array[1]}{'feature'}=$array[4];
	
}
@array=();

while(<uORF>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	#my $u = substr($array[0],0,9);
	$huORF{$array[1]}{'gene'}=$array[0];
	$huORF{$array[1]}{'mRNA'}=$array[1];
	$huORF{$array[1]}{'name'}=$array[3];
	
}
@array=();




##############Read NMD files into Hash Tables ##########
while(<U1>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	#my $u = substr($array[0],0,9);
	$hU1{$array[0]}{'mRNA'}=$array[1];
	$hU1{$array[0]}{'fc'}=$array[3];
	$hU1{$array[0]}{'pval'}=$array[4];
	
}
@array=();


while(<U3>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	$hU3{$array[0]}{'mRNA'}=$array[1];
	$hU3{$array[0]}{'fc'}=$array[3];
	$hU3{$array[0]}{'pval'}=$array[4];
	#$hVF{$array[0]}{'feat'}=$array[4];
	
	
}
@array=();


while(<UDE>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	#my $g1 = substr($array[0],0,9);
	$hUDE{$array[0]}{'gene'}=$array[0];
	$hUDE{$array[0]}{'fc'}=$array[1];
	$hUDE{$array[0]}{'pval'}=$array[3];

	
}
@array=();


while(<SDE>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	#my $g2 = substr($array[0],0,9);
	$hSDE{$array[0]}{'gene'}=$array[0];
	$hSDE{$array[0]}{'fc'}=$array[1];
	$hSDE{$array[0]}{'pval'}=$array[3];

	
}
@array=();



while(<SMGDE>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	#my $g2 = substr($array[0],0,9);
	$hSMGDE{$array[0]}{'fc'}=$array[2];
	$hSMGDE{$array[0]}{'pval'}=$array[3];

	
}
@array=();


while(<UD>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	$hUPFDE{$array[0]}{'fc'}=$array[1];

	
}
@array=();

################## Read Cap coordinates File ##########
while(<CAPS>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	$hCaps{$array[8]}{'mRNA'}=$array[11];
	$hCaps{$array[8]}{'Col_cpm'}=$array[9];
	$hCaps{$array[8]}{'xrn4_cpm'}=$array[10];

	
}
@array=();



############Read file into Hash tables ################
while(<LX>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	
	@arry=split(/\t/, $_);
	
	$hLX{$arry[0]}{'name'} = $arry[2];						
	$hLX{$arry[0]}{'col_rpkm'} = $arry[7];
	$hLX{$arry[0]}{'xrn4_rpkm'} = $arry[8];
	$hLX{$arry[0]}{'log2fc'} = $arry[9];
	$hLX{$arry[0]}{'qval'} = $arry[12];
	
}
@arry = ();

while(<LD>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	
	@arry=split(/\t/, $_);
	
	$hLD{$arry[0]}{'name'} = $arry[2];						
	$hLD{$arry[0]}{'col_rpkm'} = $arry[7];
	$hLD{$arry[0]}{'dbl_rpkm'} = $arry[8];
	$hLD{$arry[0]}{'log2fc'} = $arry[9];
	$hLD{$arry[0]}{'qval'} = $arry[12];
	
}
@arry = ();

while(<LDX>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	
	@arry=split(/\t/, $_);
	
	$hLDX{$arry[0]}{'name'} = $arry[2];						
	$hLDX{$arry[0]}{'xrn4_rpkm'} = $arry[7];
	$hLDX{$arry[0]}{'dbl_rpkm'} = $arry[8];
	$hLDX{$arry[0]}{'log2fc'} = $arry[9];
	$hLDX{$arry[0]}{'qval'} = $arry[12];
	
}
@arry = ();


while (<SEQUENCES>) {
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	
	#if($_=~ m/^(mRNA\w+)/) {
		@arry=split(/\t/, $_);	#mRNA::1:10014046-10014097(-)	CAAAACCTTCCTTTCCTAATTGGTATCTATCTTTAAAAACATACTTGAAAA
		
		my @temp = split(/::/, $arry[0]);
		my $line = $temp[1];
		my $l = length($line);
		
		my $pos = substr($line,0, $l - 3);
		#my $pos = substr($line,0, $l);
		chomp($pos);
		
		
		my $sequence = $arry[1];
		$sequence =~ tr/T/U/; #Convert T to U
		
		my $up = substr($sequence, 0, 20);
		my $down = substr($sequence, 20, 20);
		#$hSeq{$pos}{'seq'} = $sequence;
		#print "$line\t$pos\t$sequence\n";
		$hSeq{$pos}{'seq'} = lc($up) . uc($down);
		
	#}#End if	
}#End While
@arry=();	
	
while (<TABFILE>) {
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@arry=split(/\t/, $_);
	
	#Coordinates of the mRNA sequences
	my $start = $arry[3] - 1;
	
	my $locus = $arry[0].":".$start."-".$arry[4];
	#my $locus = $arry[0].":".$start."-".$arry[4] . "($arry[6])";
	chomp($locus);
	
		$hTargets{$locus}{'type'} = $arry[2];						
		$hTargets{$locus}{'start'} = $arry[1];
		$hTargets{$locus}{'ori'} = $arry[6];
		
		#Coordinates of 5'P Peak
		$hTargets{$locus}{'peak'} = $arry[8];
		
		$hTargets{$locus}{'log2_CX1'} = $arry[12]; #Log2 FC xrn4 Vs Col0 Rep1 and Rep2
		$hTargets{$locus}{'log2_CX2'} = $arry[17];
						
		$hTargets{$locus}{'log2_MX1'} = $arry[13]; #Log2 FC xrn4 Vs dbl Rep1 and Rep2
		$hTargets{$locus}{'log2_MX2'} = $arry[18];
		
		
		$hTargets{$locus}{'Col_cpm1'} = $arry[9];
		$hTargets{$locus}{'Col_cpm2'} = $arry[14];
		$hTargets{$locus}{'xrn4_cpm1'} = $arry[10];
		$hTargets{$locus}{'xrn4_cpm2'} = $arry[15];
		$hTargets{$locus}{'dbl_cpm1'} = $arry[11];
		$hTargets{$locus}{'dbl_cpm2'} = $arry[16];
		
		
		#Coordinates of the mRNA
		$hTargets{$locus}{'mRNA'} = $arry[19];		
		$hTargets{$locus}{'coord'} = $arry[20];
		$hTargets{$locus}{'Desc'} = $arry[21];
		
		
		

} #End while

my $header1 = "mRNA\tBiotype\tPeak_Pos\tPeak_Seq\tPeak_Region\tdne1_Dependent_Site\tLocation\tMajor_Site_in_xrn4\t";
my $header0 = "Col0_cpm1\txrn4_cpm1\tdne1xrn4_cpm1\tLog2FC_xrn4_v_Col0_Rep1\tLog2FC_xrn4_v_dne1xrn4_Rep1\t";
my $header2 = "Col0_cpm2\txrn4_cpm2\tdne1xrn4_cpm2\tLog2FC_xrn4_v_Col0_Rep2\tLog2FC_xrn4_v_dne1xrn4_Rep2\tmRNA\tLocus\tAnnotation\tName\t";
my $header3 = "Col0_rpkm\txrn4_rpkm\tdne1xrn4_rpkm\tLog2FC_xrn4_v_Col0\tq-value\tLog2FC_dne1xrn4_v_Col0\tq-value\tLog2FC_dne1xrn4_v_xrn4\tq-value\t";
my $header4 = "Cap_Pos\tCol0_cpm\txrn4_cpm\tDecapped\t";
my $header5 = "Gene\tupf1pad4_log2FC\tupf1pad4_reg\tsmg7pad4_log2FC\tsmg7pad4_reg\tupf1_FC\tupf3_FC\tupf1upf3_FC\tsmg7_pad4_log2FC\tsmg7_pad4_qval\n";
my $destination = "$header1" . "$header0" . "$header2" . "$header3" . "$header4" . "$header5";

print OUT $destination;

my $posi;
my $peakpos;
my $transcript;
my $gid;
my ($cap1, $cap2, $reg) = "";

foreach my $key (sort keys %hTargets){
	$posi = $key;
	$peakpos = $hTargets{$key}{'peak'};
	$transcript = $hTargets{$key}{'mRNA'};
	chomp($transcript);
	$gid = substr($transcript,0,9);
	#my $cappos = $hCaps{$key}{'mRNA'};
	$cap1 = $hCaps{$peakpos}{'Col_cpm'};
	$cap2 = $hCaps{$peakpos}{'xrn4_cpm'};
	
	
	if($cap1 >= 1 && $cap2 >= 1) {
			$reg = "Decapped";
	} else {
		$reg = "";
	}
	
	my $transcript = $hTargets{$key}{'mRNA'};
	my $site = $hTargets{$key}{'peak'};
	
    print OUT "$hTargets{$key}{'mRNA'}\t$hTargets{$key}{'type'}\t$hTargets{$key}{'peak'}\t$hSeq{$posi}{'seq'}\t$key\t";
	print OUT "$hSites{$site}{'peak_type'}\t$hSites{$site}{'feature'}\t$hSites{$site}{'major_peak'}\t";
	print OUT "$hTargets{$key}{'Col_cpm1'}\t$hTargets{$key}{'xrn4_cpm1'}\t$hTargets{$key}{'dbl_cpm1'}\t$hTargets{$key}{'log2_CX1'}\t$hTargets{$key}{'log2_MX1'}\t";
    print OUT "$hTargets{$key}{'Col_cpm2'}\t$hTargets{$key}{'xrn4_cpm2'}\t$hTargets{$key}{'dbl_cpm2'}\t$hTargets{$key}{'log2_CX2'}\t$hTargets{$key}{'log2_MX2'}\t";
	print OUT "$hTargets{$key}{'mRNA'}\t$hTargets{$key}{'coord'}\t$hTargets{$key}{'Desc'}\t$hLX{$transcript}{'name'}\t";
	print OUT "$hLX{$transcript}{'col_rpkm'}\t$hLX{$transcript}{'xrn4_rpkm'}\t$hLDX{$transcript}{'dbl_rpkm'}\t$hLX{$transcript}{'log2fc'}\t$hLX{$transcript}{'qval'}\t$hLD{$transcript}{'log2fc'}\t$hLD{$transcript}{'qval'}\t";
	print OUT "$hLDX{$transcript}{'log2fc'}\t$hLDX{$transcript}{'qval'}\t";
	print OUT "$peakpos\t$cap1\t$cap2\t$reg\t";
	print OUT "$gid\t$hUDE{$gid}{'fc'}\t$hUDE{$gid}{'pval'}\t$hSDE{$gid}{'fc'}\t$hSDE{$gid}{'pval'}\t$hU1{$gid}{'fc'}\t$hU3{$gid}{'fc'}\t$hUPFDE{$gid}{'fc'}\t$hSMGDE{$gid}{'fc'}\t$hSMGDE{$gid}{'pval'}\t$huORF{$transcript}{'mRNA'}\n";
	

    
}
close(OUT);

###########################################################