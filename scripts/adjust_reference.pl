#!/usr/bin/perl
use strict;

my $vcf = $ARGV[0]; #VCF containing SNPs to be tested
my $id = $ARGV[1]; #ID in vcf

open (OUT, ">map_over.txt") || die "Unable to open file to write to: $!\n";

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
    my @array = split;
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
		print OUT "$array[0]\t$array[1]\t$array[3]\t$array[4]\t$array[9]\t$f\t$m\n";
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
close (OUT);
    
