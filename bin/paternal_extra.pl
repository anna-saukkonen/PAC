#!/usr/bin/perl

open (F,"paternal_tags_UM.txt") or die $!; my %tag=(); while (<F>) { chomp; $tag{$}=0; } close (F); open (N,"paternal_tags_UM.RSEM.txt") or die $!; while (<N>) { chomp; if (!(exists $tag{$})) { print "$_\n"; }} close (N);
