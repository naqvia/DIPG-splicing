#!/usr/bin/perl

my @splice_ids;
my (%splice_perc_wt,%splice_perc_mt);

while(<>)
{
  chomp;
  next if($_=~/SpliceID/);
  my ($id,$freq,$prev,$type) = split ",";
  push @splice_ids,$id;

  if($type=~/WT/)
  {
    $splice_perc_wt{$id} = $prev;
  }
  else{
    $splice_perc_mt{$id} = $prev;
  }
}

my @splice_ids_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_ids };

foreach my $splice_id(@splice_ids_uniq)
{
  print $splice_id,",";
  if($splice_perc_wt{$splice_id}) {
    print $splice_perc_wt{$splice_id};
  }
  else{
    print "0";
  }
  print ",";

  if($splice_perc_mt{$splice_id}) {
    print $splice_perc_mt{$splice_id};
  }
  else{
    print "0";
  }
  print "\n";
}
__DATA__
SpliceID,Freq,Prev,Type
AACS_chr12:125136661-125136864,7,0.636364,WT
AAK1_chr2:69476879-69476990,3,0.272727,WT
AAK1_chr2:69495984-69496080,8,0.727273,WT
AAMDC_chr11:77841018-77841060,8,0.727273,WT
AASS_chr7:122082787-122082891,3,0.272727,WT
ABCA11P_chr4:437445-437511,4,0.363636,WT
ABCA2_chr9:137023837-137023840,6,0.545455,WT
ABCA5_chr17:69251746-69251866,2,0.181818,WT
ABCA8_chr17:68949311-68949472,1,0.0909091,WT
