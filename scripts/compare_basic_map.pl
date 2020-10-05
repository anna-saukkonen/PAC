use strict;

#name: compare_basic_map.pl
#description: A program to obtain reads and pairs for heterozygote positions in RNA data from two parental genome mappings

my $snp_file = $ARGV[0]; #snp to test file
my $bam = $ARGV[1]; #final mother bam
my $id = $ARGV[2]; #VCF ID
my $res1 = $ARGV[3]; #baq res
my $res2 = $ARGV[4]; #non-baq res

my %ref = ();
my %alt = ();

open (OUTM, ">extract.bed") || die "Unable to open mother file to write to: $!\n";


my %hets = (); my @het_array = ();
my @head = ();
open (EXOME, "$snp_file") || die "Unable to open snp file to read to: $!\n";
while (<EXOME>) {  #get edit positions and edit types
  if ($_ =~ /^#CHROM/) {
    @head = split;
  } 
  if (!($_ =~ /^#CHROM/)) {
    my @array = split;
    for (my $i=0;$i<@array;$i++) {
      if ($head[$i] eq "$id") {
	if (($array[$i] =~ /0[\/\|]1/)||($array[$i] =~ /1[\/\|]0/)) {
	  my $chr = $array[0]; $chr =~ s/chr//; $chr = "chr".$chr; #Make sure that vcf contains "chr";
	  my $tag = $chr."_".$array[1];
	  if (($array[3] =~ /^[AGTC]$/)&&($array[4] =~ /^[AGTC]$/)) { #Single bp SNPs only
	    if (!($array[4] =~ /,/)) {  #just look at single allele sites (may need to improve this later)
	      $hets{$tag} = $array[3]."_".$array[4]; #store het positions with base change
	      push @het_array, $tag;
	      my $before=$array[1]-5; my $after=$array[1]+5; 
	      print OUTM "$chr\t$before\t$after\n";
	      $ref{$tag}=$array[3];
	      $alt{$tag}=$array[4];
	    }
	  }
	}
      }
    }
  }
}
close (EXOME);
close (OUTM);

print "SNPs\n";

system ("samtools mpileup -x -d 10000000 -l extract.bed $bam > pileup_baq.pu");
system ("samtools mpileup -x -B -d 10000000 -l extract.bed $bam > pileup.pu");

##Now: run through pileup and store in hash the comparible positions for mother and father
##Then, run through hash and convert refs to dots - then pump into subroutine and produce output.

my %combined_hash = ();

open (PILEUP1, "pileup_baq.pu") || die "Unable to open real pileup file to read: $!\n";
while (<PILEUP1>) {
  my @array = split;
  my $tag = $array[0]."_".$array[1];
  if (exists $hets{$tag}) {
    $combined_hash{$tag} .= $array[4];
  }
}
close (PILEUP1);

open (RES, ">$res1") || die "Unable to open results file to write to: $!\n";
print RES "#Chr\tPos\tRefAl\tAltAl\tMapRef\tMapAlt\tMapRatio\tMapCov\n";

foreach (@het_array) {
  my $seq = $combined_hash{$_};
  #print "$_\t$seq\n";
  $seq =~ s/$ref{$_}/\./gi;
  my ($ref,$alt,$prop,$cov) = &countpos($seq,$seq,$alt{$_});
  my @a=split(/_/,$_);
  print RES "$a[0]\t$a[1]\t$ref{$_}\t$alt{$_}\t$ref\t$alt\t$prop\t$cov\n" if ($cov>0);
}

close (RES);

%combined_hash = ();

open (PILEUP3, "pileup.pu") || die "Unable to open real pileup file to read: $!\n";
while (<PILEUP3>) {
  my @array = split;
  my $tag = $array[0]."_".$array[1];
  if (exists $hets{$tag}) {
    $combined_hash{$tag} .= $array[4];
  }
}

close (PILEUP3);

open (RES1, ">$res2") || die "Unable to open results file to write to: $!\n";
print RES1 "#Chr\tPos\tRefAl\tAltAl\tMapRef\tMapAlt\tMapRatio\tMapCov\n";

foreach (@het_array) {
  my $seq = $combined_hash{$_};
  #print "$_\t$seq\n";
  $seq =~ s/$ref{$_}/\./gi;
  my ($ref,$alt,$prop,$cov) = &countpos($seq,$seq,$alt{$_});
  my @a=split(/_/,$_);
  print RES1 "$a[0]\t$a[1]\t$ref{$_}\t$alt{$_}\t$ref\t$alt\t$prop\t$cov\n" if ($cov>0);
}

close (RES1);

#system ("rm mother.bed");
#system ("rm father.bed");
#system ("pat_pileup_baq.pu");
#system ("pat_pileup.pu");
#system ("mat_pileup_baq.pu");
#system ("mat_pileup.pu");

###subroutines

sub countpos ($$$$) {
  my ($seq,$score,$alt,$refb) = @_;
  $seq =~ s/\$//g; #end of read tag - not associated with any quality score
  $seq =~ s/\^.//g; #start of read tag - always followed by a mapping quality for the read, this is different from the base quality in the final column
  my @s = split (//,$seq); my $r = 0; my $c = 0;
  for (my $i=0;$i<@s;$i++) {
    if ($s[$i] =~ /[\.\,$alt]/i) {
      $c++;
      if ($s[$i] =~ /[\.\,]/i) {
	$r++;
      }
    }
    if ($s[$i] =~ /[\d]/) {
      $i+=$s[$i];
    }
  }
  my $p = 0; my $a = 0;
  if ($c>0) {
    $p = $r/$c;
    $a = $c-$r;
  }
  return($r,$a,$p,$c);
}

