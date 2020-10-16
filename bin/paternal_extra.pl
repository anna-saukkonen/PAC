#!/usr/bin/perl

open (F,"paternal_tags_UM.txt"); my %tag=(); while (<F>) { chomp; $tag{$_}=0; } close (F); open (N,"paternal_tags_UM.RSEM.txt"); while (<N>) { chomp; if (!(exists $tag{$_})) { print "$_\n"; }} close (N);
