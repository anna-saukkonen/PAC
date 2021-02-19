#!/usr/bin/perl
use strict;

my $bed = $ARGV[0]; #original gene bed file
my $mother_ae = $ARGV[1]; #mother ae file
my $father_ae = $ARGV[2]; #father ae file
my $id = $ARGV[3]; #ID for individual

my $out = $id."_gene_level_ae.txt";

open (OUT, ">$out") || die "Unable to open file to write to: $!\n";

my %gene_start=(); my %gene_end=(); my %gene_chr=();
open (BED, "$bed") || die "Unable to open bed file: $!\n";
while (<BED>) {
  my @array=split;
  $gene_chr{$array[3]}=$array[0];
  $gene_start{$array[3]}=$array[1];
  $gene_end{$array[3]}=$array[2];
}
close (BED);

open (MOTH, "$mother_ae") || die "Unable to open mother ae file: $!\n";
my %moth_a=(); my %moth_b=();
while (<MOTH>) {
  my @array=split;
  if (!($_=~/^contig/)) {
    $moth_a{$array[3]}=$array[4];
    $moth_b{$array[3]}=$array[5];
  }
}
close (MOTH);

print OUT "contig\tstart\tstop\tname\taCount\tbCount\ttotalCount\n";

open (FATH, "$father_ae") || die "Unable to open father ae file: $!\n";
while (<FATH>) {
  my @array=split;
  if (!($_=~/^contig/)) {
    my $a=$moth_a{$array[3]}+$array[4];
    my $b=$moth_b{$array[3]}+$array[5];
    my $tot=$a+$b;
    print OUT "$gene_chr{$array[3]}\t$gene_start{$array[3]}\t$gene_end{$array[3]}\t$array[3]\t$a\t$b\t$tot\n";
  }
}
close (FATH);
close (OUT);
