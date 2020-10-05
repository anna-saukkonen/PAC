use strict;

###Created By: Alan Hodgkinson
###Program to select best macth for a read across two BAM files. Score: +1 for correct match, -2 for indel.

my $maternal = $ARGV[0]; #Maternal sam file;
my $paternal = $ARGV[1]; #Paternal sam file;

open (MAT, "$maternal") || die "Unable to open maternal file to read: $!\n";

my @maternal = ();
my $mark = 0; my $check = 0;;

my $count=0; my $score=0; my $tag; my $old_tag;
while (<MAT>) {
  $mark++;
  if ((int($mark/100000))>$check) { 
    print "$check\n";
    $check = int($mark/100000);
  }
  chomp;
  my @array = split; 
  $old_tag=$tag;
  $tag = $array[0];
  my $match;
  foreach (@array) {
    $match = &ParseMDscore($_) if ($_ =~ /MD/);
  }
  my ($indel,$cmatch) = &ParseCigar($array[5]);
  $score+=$match; $score-=(2*$indel);
  $count++;
  if ($count==2) {
    my $final = $tag."splitting".$score;
    push @maternal, $final;
    $count=0; $score=0;
    print "Tags not paired\n" if ($tag ne $old_tag);
  }
}
    
close (MAT);

open (PAT, "$paternal") || die "Unable to open paternal file to read: $!\n";
open (OUTMAT, ">maternal_wins.txt") || die "Unable to open maternal file to write to: $!\n";
open (OUTPAT, ">paternal_wins.txt") || die "Unable to open paternal file to write to: $!\n";


my $point=0;
$count=0; $score=0; $tag; $old_tag;
while (<PAT>) {
  chomp;
  my @array = split;
  $old_tag=$tag;
  $tag = $array[0];
  my $match;
  foreach (@array) {
    $match = &ParseMDscore($_) if ($_ =~ /MD/);
  }
  my ($indel,$cmatch) = &ParseCigar($array[5]);
  $score+=$match; $score-=(2*$indel);
  $count++;
  if ($count==2) {
    print "Tags not paired\n" if ($tag ne $old_tag);
    my $mfinal = $maternal[$point];
    my @marray=split(/splitting/,$mfinal);
    my $mtag = $marray[0];
    my $mscore = $marray[1];
    warn "Tags don't match: $tag\t$mtag\n" if ($tag ne $mtag);
    print OUTPAT "$tag\n" if ($score>$mscore);
    print OUTMAT "$tag\n" if ($score<$mscore);
    if ($score==$mscore) {
      my $r = int(rand(2));
      print OUTPAT "$tag\n" if ($r==0);
      print OUTMAT "$tag\n" if ($r==1);
    }
    $count=0; $score=0; $point++;
  }
}
    
close (PAT);

#############
#Subroutines#
#############

sub ParseMDscore ($) {
  my $field = @_[0]; #get MD tag
  my @a = split (/:/,$field);
  my $string = $a[2]; my @scan = split (//,$string); #collect info section and split
  my $sum=0; my $store;
  foreach (@scan) { 
    if ($_ =~ /\d/) { #if a digit
      $store .= $_; #add to store string
    }
    if ($_ =~ /\D/) { #if non-digit
      $sum+=$store; #add to sum (thus collecting sum of correclty matching nucs)
      $store = "";
    }
  }
  $sum+=$store; #at end, add any remaining stored numbers to sum
  return($sum);
}

sub ParseCigar ($) {
  my $field = @_[0];
  my @scan = split (//,$field); #split cigar field
  my $sum=0; my $sum_m=0; my $store;
  foreach (@scan) {
    if ($_ =~ /\d/) { #if digit, store
      $store .= $_;
    }
    if ($_ =~ /\D/) { #if non-digit
      if ($_ =~ /[ID]/) { #if indel, count
	$sum+=$store;
      }
      if ($_ =~ /M/) { #if match, count
	$sum_m+=$store;
      }
      $store = ""; #clear store
    }
  }
  return($sum,$sum_m);
}
