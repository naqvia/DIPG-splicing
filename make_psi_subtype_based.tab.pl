#!/usr/bin/perl

#
#
#

my @files  = <majiq_voila_tsv/*deltapsi.voila.tsv>;
my (%lsv_gene,%psi_lsv);
my  (@samples,@lsvs);
my %H3K28_bs;

foreach my $file (@files)
{
  my $sample = $file;
  if($file=~/(BS\_[\w+\d+]+)\.d/)
  {
    #print "sample: ",$1,"\t";
    $sample = $1;
  }

  push @samples,$sample;

  open(FIL,$file) ||  die("Cannot Open File");
  while(<FIL>)
  {
    my @cols = split "\t";
    my $gene= $cols[0];
    my $lsv = $cols[2];
    my $psi = $cols[7];
    my (@psis) = split/\;/,$psi;
    $first_psi  = $psis[0];
    push @lsvs,$lsv;
    $lsv_gene{$lsv} = $gene;
    #print "lsv:".$lsv,"\tsample:",$sample,"\tpsi:",$psi,"_".$first_psi,"\n";
    $psi_lsv{$lsv}{$sample} = $psis[0];

  }
}

my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @lsvs_uniq = do { my %seen; grep { !$seen{$_}++ } @lsvs };
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


print "lsv";
foreach my $sample(@samples_uniq)
{

  my $sample_to_check = $sample;
  $sample_to_check=~s/BS_tumor_//;

  unless ($H3K28_bs{$sample_to_check})
  {
    print ",";
    print $sample."_Non_H3K28";
  }
}
print ",";

foreach my $sample(@samples_uniq)
{
  my $sample_to_check2 = $sample;
  $sample_to_check2=~s/BS_tumor_//;
  if  ($H3K28_bs{$sample_to_check2})
  {
    print ",";
    print $sample."_H3K28";
  }
}
print "\n";

foreach my $lsv(@lsvs_uniq)
{
  print $lsv_gene{$lsv}."_".$lsv;
  foreach my $sample(@samples_uniq)
  {
    my $sample_to_check = $sample;
    $sample_to_check=~s/BS_tumor_//;
    unless  ($H3K28_bs{$sample_to_check})
    {
      print ",";
      #print "sample:",$sample.":".$lsv.":";
      if($psi_lsv{$lsv}{$sample} > 0){
        print $psi_lsv{$lsv}{$sample};
      }
      else{
        print "0";
      }
    }
  }
  foreach my $sample(@samples_uniq)
  {
    my $sample_to_check = $sample;
    $sample_to_check=~s/BS_tumor_//;
    if  ($H3K28_bs{$sample_to_check})
    {
      print ",";
      #print "sample:",$sample.":".$lsv.":";
      if($psi_lsv{$lsv}{$sample}){
        print $psi_lsv{$lsv}{$sample};
      }
      else{
        print "0";
      }
    }
  }
  print "\n";

}
__DATA__
TOM1L2	ENSG00000175662.17:s:17861476-17861551	0.845
