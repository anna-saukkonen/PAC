#!/usr/bin/perl
use strict;

my $parent = $ARGV[0]; #paternal or maternal

my $file1 = $parent."_tags_UM.txt";
my $file2 = $parent."_tags_UM.RSEM.txt";
my $out = "extra.rsem.".$parent.".txt";
  
open (F,"$file1");
my %tag=();

while (<F>) {
  chomp; $tag{$_}=0;
}
close (F);

open (N,"$file2");
open (OUT, ">$out") || die "Unable to open file to write to: $!\n";

while (<N>) {
  chomp;
  if (!(exists $tag{$_})) {
    print OUT "$_\n";
  }
}

close (N);
close (OUT);
