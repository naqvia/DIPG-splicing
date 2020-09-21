#!/usr/bin/perl

#
#
#

my @files  = <majiq_voila_tsv/*deltapsi.voila.tsv>;
my (%lsv_gene,%psi_lsv);
my  (@samples,@lsvs);

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

print "lsv";
foreach my $sample(@samples_uniq)
{
  print "\t";
  print $sample;
}
print "\n";

foreach my $lsv(@lsvs_uniq)
{
  print $lsv_gene{$lsv}."_".$lsv;
  foreach my $sample(@samples_uniq)
  {
    print "\t";
    if($psi_lsv{$lsv}{$sample}>0){
      print $psi_lsv{$lsv}{$sample};
    }
    else{
      print "0";
    }
  }
  print "\n";
}

__DATA__
TOM1L2	ENSG00000175662.17:s:17861476-17861551	0.845
