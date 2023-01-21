#!/usr/bin/perl -w
#use strict;

my $sites = "TS1_Leaf.txt"; #MaxSeq and Major internal sites finalized
my $cds = "TAIR10_v49_CDS_Tabsep.txt"; #Tabsep CDS sequence file
my $threeutr = "TAIR10_3UTR_Tabsep.txt"; #Tabsep 3UTR sequence file
my $fiveutr = "TAIR10_5UTR_Tabsep.txt"; #Tabsep 5UTR sequence file

my $output = "TS1_Leaf_MetaGene_Profile.txt";

open(CDS, $cds) || die "Can't open $cds";
open(UTR3, $threeutr) || die "Can't open $threeutr";
open(UTR5, $fiveutr) || die "Can't open $fiveutr";

open(SITES, $sites) || die "Can't open $sites";


#######################################################

open (OUT, ">$output") || die "Can't open $output";


#Read headers
<SITES>; #read headers

my %hTargets = ();
my %hCDS = ();
my %hUTR3 = ();
my %hUTR5 = ();

my @array=();

####################### Read Final site list #############
while(<SITES>){
	chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	
	$hTargets{$array[0]}{'coord'}=$array[1];

	my $seq = substr($array[2], 20, 20);
	$seq =~ tr/U/T/;
	#$seq =~ tr/u/t/;
	
	$hTargets{$array[0]}{'peak_seq'}=$seq;
	$hTargets{$array[0]}{'feature'}=$array[3]; #Intron / Exon / 3'UTR and 5'UTR
		
	$hTargets{$array[0]}{'xrn4_cpm'}=$array[5];
	$hTargets{$array[0]}{'col_cpm'}=$array[4];
	$hTargets{$array[0]}{'dbl_cpm'}=$array[6];

}
@array=();

while(<CDS>){
chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	
	$hCDS{$array[0]}{'seq'}=$array[3];

}

@array=();

while(<UTR3>){
chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	
	$hUTR3{$array[0]}{'seq'}=$array[3];

}
@array=();

while(<UTR5>){
chomp $_;
	s/\r//; #hidden character removal
	s/\f//; #hidden character removal
	@array=split(/\t/, $_);
	
	$hUTR5{$array[0]}{'seq'}=$array[3];

}

#######################################################
print OUT "Transcript\tSequence\tFeature\tScore\tCol0_cpm\txrn4_cpm\tmfl1xrn4_cpm\n";

foreach my $key (sort keys %hTargets){

my $feature = $hTargets{$key}{'feature'};
my $sequence = $hTargets{$key}{'peak_seq'};
my $source;
my $metascore;



if ($feature eq "Exon"){
	#Find position of seq 
	$source = $hCDS{$key}{'seq'};
	$metascore = getScore($sequence, $source);
	$metascore = $metascore + 1; #Score 1.1 to 2.0
	
}elsif ($feature eq "3-UTR"){
	$source = $hUTR3{$key}{'seq'};
	$metascore = getScore($sequence, $source);
	$metascore = $metascore + 2; #Score 2.1 to 3
	print "$feature\t$sequence\t$source\n";

}elsif ($feature eq "5-UTR"){
	$source = $hUTR5{$key}{'seq'};
	$metascore = getScore($sequence, $source);
	$metascore = $metascore; #Score 0 to 1

}else {
# No feature found
	$metascore = "NA";

}

print OUT "$key\t$sequence\t$feature\t$metascore\t$hTargets{$key}{'col_cpm'}\t$hTargets{$key}{'xrn4_cpm'}\t$hTargets{$key}{'dbl_cpm'}\n";


}


sub getScore {
my $seq =  $_[0];
my $src =  $_[1];

my $elen = length($src);

if($elen > 0){
	my $maxpos = index($src, $seq) + 1;
	my $fract = 1 - (($elen - $maxpos)/$elen);
	return $fract;
	print $fract;
	print "\n";
}


} #End sub

