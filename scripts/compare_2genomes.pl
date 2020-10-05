use strict;

#name: compare_2genomes.pl
#description: A program to obtain reads and pairs for heterozygote positions in RNA data from two parental genome mappings

my $map = $ARGV[0]; #adjust ref file
my $snp_file = $ARGV[1]; #snp to test file
my $mother_bam = $ARGV[2]; #final mother bam
my $father_bam = $ARGV[3]; #final father bam
my $id = $ARGV[4]; #VCF ID
my $res1 = $ARGV[5]; #baq res
my $res2 = $ARGV[6]; #non-baq res

my %mother_swap = ();
my %father_swap = ();
my %ref = ();
my %alt = ();
open (SWAP, "$map") || die "Unable to open file to read: $!\n";
open (OUTM, ">mother.bed") || die "Unable to open mother file to write to: $!\n";
open (OUTF, ">father.bed") || die "Unable to open mother file to write to: $!\n";

while (<SWAP>) {
  my @array = split;
  my $chr=$array[0]; $chr=~s/chr//;
  my $tag = "chr".$chr."_".$array[1];
  my $tagf = "chr".$chr."_".$array[5];
  my $tagm = "chr".$chr."_".$array[6];
  $ref{$tag} = $array[2];
  $alt{$tag} = $array[3];
  $father_swap{$tagf} = $tag;
  $mother_swap{$tagm} = $tag;
  my $fb=$array[5]-1; my $fa = $array[5]+1;
  print OUTF "chr$chr\_paternal\t$fb\t$fa\n";
  my $mb=$array[6]-1; my $ma = $array[6]+1;
  print OUTM "chr$chr\_maternal\t$mb\t$ma\n";
}

close (SWAP);
close (OUTM);
close (OUTF);

print "SWAP\n";

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
	    }
	  }
	}
      }
    }
  }
}
close (EXOME);

print "SNPs\n";

system ("samtools mpileup -x -d 10000000 -l father.bed $father_bam > pat_pileup_baq.pu");
system ("samtools mpileup -x -B -d 10000000 -l father.bed $father_bam > pat_pileup.pu");
system ("samtools mpileup -x -d 10000000 -l mother.bed $mother_bam > mat_pileup_baq.pu");
system ("samtools mpileup -x -B -d 10000000 -l mother.bed $mother_bam > mat_pileup.pu");

##Now: run through pileup and store in hash the comparible positions for mother and father
##Then, run through hash and convert refs to dots - then pump into subroutine and produce output.

my %combined_hash = ();

open (PILEUP1, "pat_pileup_baq.pu") || die "Unable to open real pileup file to read: $!\n";
while (<PILEUP1>) {
  my @array = split;
  my $tag = $array[0]."_".$array[1];
  $tag =~ s/_paternal//;
  my $new_tag = $father_swap{$tag};
  if (exists $hets{$new_tag}) {
    $combined_hash{$new_tag} .= $array[4];
  }
}
close (PILEUP1);

open (PILEUP2, "mat_pileup_baq.pu") || die "Unable to open real pileup file to read: $!\n";
while (<PILEUP2>) {
  my @array = split;
  my $tag = $array[0]."_".$array[1];
  $tag =~ s/_maternal//;
  my $new_tag = $mother_swap{$tag};
  if (exists $hets{$new_tag}) {
    $combined_hash{$new_tag} .= $array[4];
  }
}
close (PILEUP2);

open (RES, ">$res1") || die "Unable to open results file to write to: $!\n";
print RES "#Chr\tPos\tRefAl\tAltAl\tMapRef\tMapAlt\tMapRatio\tMapCov\n";

foreach (@het_array) {
  my $seq = $combined_hash{$_};
  if ($seq=~/[agtc\.\,]/i) {
    #print "$_\t$seq\n";
    $seq =~ s/$ref{$_}/\./gi;
    my ($ref,$alt,$prop,$cov) = &countpos($seq,$seq,$alt{$_});
    my @a=split(/_/,$_);
    print RES "$a[0]\t$a[1]\t$ref{$_}\t$alt{$_}\t$ref\t$alt\t$prop\t$cov\n";
  }
}

close (RES);

%combined_hash = ();

open (PILEUP3, "pat_pileup.pu") || die "Unable to open real pileup file to read: $!\n";
while (<PILEUP3>) {
  my @array = split;
  my $tag = $array[0]."_".$array[1];
  $tag =~ s/_paternal//;
  my $new_tag = $father_swap{$tag};
  if (exists $hets{$new_tag}) {
    $combined_hash{$new_tag} .= $array[4];
  }
}

close (PILEUP3);

open (PILEUP4, "mat_pileup.pu") || die "Unable to open real pileup file to read: $!\n";
while (<PILEUP4>) {
  my @array = split;
  my $tag = $array[0]."_".$array[1];
  $tag =~ s/_maternal//;
  my $new_tag = $mother_swap{$tag};
  if (exists $hets{$new_tag}) {
    $combined_hash{$new_tag} .= $array[4];
  }
}
close (PILEUP4);

open (RES1, ">$res2") || die "Unable to open results file to write to: $!\n";
print RES1 "#Chr\tPos\tRefAl\tAltAl\tMapRef\tMapAlt\tMapRatio\tMapCov\n";

foreach (@het_array) {
  my $seq = $combined_hash{$_};
  if ($seq=~/[agtc\.\,]/i) {
    #print "$_\t$seq\n";
    $seq =~ s/$ref{$_}/\./gi;
    my ($ref,$alt,$prop,$cov) = &countpos($seq,$seq,$alt{$_});
    my @a=split(/_/,$_);
    print RES1 "$a[0]\t$a[1]\t$ref{$_}\t$alt{$_}\t$ref\t$alt\t$prop\t$cov\n";
  }
}

close (RES1);

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

