#!/usr/bin/perl

#
# make_tab_for_global_splicing.pl
#

my @tsv = <majiq_voila_tsv/dipg_deltapsinormal_control-BS_tumor_BS_*tsv>;
my %lsv_counts;
my %H3K28_samples;
my (%lsv_counts_mutant, %lsv_counts_other);
my %lsv_counts_total;

open(FIL,"dipg_samples_w_molecsubtype.txt");
while(<FIL>)
{
  chomp;
  my ($bs, $pt, $subtype) = split "\t";
  if($subtype=~/H3\sK28/){
    $H3K28_bs{$bs} = $bs;
    #$H3K28_bs{"H3K28"} = $bs;
  }
}

foreach my $tsv (@tsv)
{
  my $sample = "";
  if($tsv=~/(BS\_[\w+\d+]+)\.d/)
  {
    #print "sample: ",$1,"\t";
    $sample = $1;
    $sample=~s/BS_tumor_//;
  }

  open(FIL, $tsv) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;
    my @cols = split "\t";
    my $gene = $cols[0];
    my $lsv = $cols[2];
    #print "sample:",$sample,"\n";
    if($H3K28_bs{$sample})
    {
      $lsv_counts_mutant{$lsv}++;
    }
    else{
      $lsv_counts_other{$lsv}++;
    }
    $lsv_counts_total{$lsv}++;
  }
}

foreach my $lsv_mut(keys %lsv_counts_mutant)
{
  print "H3K28,".$lsv_mut,",",$lsv_counts_mutant{$lsv_mut},"\n";
}

foreach my $lsv_other(keys %lsv_counts_other)
{
  print "Non_H3K28,".$lsv_other,",",$lsv_counts_other{$lsv_other},"\n";
}
foreach my $lsv(keys %lsv_counts_total)
{
  print "total,".$lsv,",",$lsv_counts_total{$lsv},"\n";
}
