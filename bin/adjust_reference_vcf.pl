#!/usr/bin/perl
use strict;

my $vcf = $ARGV[0]; #VCF containing SNPs to be tested
my $id = $ARGV[1]; #ID in vcf
my @vcf=split(/\//,$vcf);
my $mother_vcf=pop(@vcf);
my $father_vcf=$mother_vcf;
$mother_vcf=~s/vcf//; $mother_vcf.="mother.vcf";
$father_vcf=~s/vcf//; $father_vcf.="father.vcf";

open (OUTM, ">$mother_vcf") || die "Unable to open file to write to: $!\n";
open (OUTF, ">$father_vcf") || die "Unable to open file to write to: $!\n";

for (1..24) {
  my $chr = $_;
  $chr="X" if ($chr==23);
  $chr="Y" if ($chr==24);
  my $file = "chr".$chr."_".$id.".map";
  my @reference = (); my @father = (); my @mother = ();
  open (MAP, "$file") || die "Unable to open map file $file to read: $!\n";
  my $m=0; my $f=0;
  while (<MAP>) {
    my @array = split;
    if ($array[0]>0) {
      push @reference, $array[0];
      if ($array[1]>0) {
	$f = $array[1]-$array[0];
      }
      if ($array[2]>0) {
	$m = $array[2]-$array[0];
      }
      push @father, $f;
      push @mother, $m;
    }
  }
  close (MAP);
  
  push @reference, 1000000000;
  
  open (SNPS, "$vcf") || die "Unable to open SNP file to read: $!\n";
  my $point=0; my @head = ();
  while (<SNPS>) {
    chomp;
    my @array = split;
    if (($_ =~ /^#/)&&($chr==1)) {
      print OUTF "$_\n";
      print OUTM "$_\n";
    }
    if ($_ =~ /^#CHROM/) {
      @head = split;
    }
    for (my $i=0;$i<@array;$i++) {
      if ($head[$i] eq "$id") {
	if (($array[$i] =~ /0[\/\|]1/)||($array[$i] =~ /1[\/\|]0/)) {
	  if ($array[0] =~ /^(chr)?$chr$/) {
	    if ((length($array[3])==1)&&(length($array[4])==1)) {
	    redoloop:
	      if (($array[1]>$reference[$point])&&($array[1]<$reference[$point+1])) {
		my $m=$array[1]+$mother[$point];
		my $f=$array[1]+$father[$point];
		print OUTF "$array[0]_paternal\t$f";
		print OUTM "$array[0]_maternal\t$m";
		for (my $j=2;$j<@array;$j++) {
		  print OUTF "\t$array[$j]";
		  print OUTM "\t$array[$j]";
		}
		print OUTF "\n";
		print OUTM "\n";
	      }
	      if (($array[1]>$reference[$point+1])) {
		$point++;
		goto redoloop;
	      }
	    }
	  }
	}
      }
    }
  }
  close (SNPS);
}
close (OUTF);
close (OUTM);

