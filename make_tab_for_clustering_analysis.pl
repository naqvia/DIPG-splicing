#!/usr/bin/perl

##
# make_tab_for_clustering_analysis.pl
#

my $dir_file  = $ARGV[0];
my $flag_progr = $ARGV[1];
my %spliced_event;
my %count_spliced_event;
my (@sample_prog, @sample_nonprog);

if($flag_progr eq 1 )
{
  open(FIL,"/Users/naqvia/Desktop/DIPG/dipg_samples_w_prog.txt");
  while(<FIL>)
  {
    chomp;
    my ($sample,$descr) = split "\t";
    if($descr = ~/Prog/)
    {
      push @sample_prog,$sample;

    }
    else{
      push @sample_nonprog, $sample;
    }
  }
}
open(FIL,$dir_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my $file = $_;
  my $sample = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $sample = $1;
  }
  open(SE_FIL,$file);
  while(<SE_FIL>)
  {
    chomp;
    if($_=~/^\d+/)
    {
      #print "*".$_,"\n";

      my @cols = split "\t";
      my $gene = $cols[2];
      my $chr  = $cols[3];
      my $str  =   $cols[4];

      my $exon_start    = $cols[5];
      my $exon_end      = $cols[6];
      my $tumor_IJC     = $cols[14];
      my $tumor_SJC     = $cols[15];
      my $inc_from_len  = $cols[16];
      my $skip_from_len = $cols[17];
      my $pval = $cols[18];
      my $thr_diff = $cols[-1];
      #next unless (abs($thr_diff) >= .10);
      next unless ($pval<=0.05);
      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);


      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end;

      $splice_id=~s/\"//g;
      push @samples, $sample;
      push @splice_id, $splice_id;
      $samples_filter{$sample} = $sample;

      #print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";
      #print $splice_id,"\t",$thr_diff,"\n";

      $spliced_event{$splice_id}{$sample} =  $thr_diff;
      push @{$count_spliced_event{$splice_id}},$thr_diff;

    }
  }
  close(SE_FIL);
}

my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_id };




print "SpliceID";
foreach my $sample(@samples_uniq)
{
  print ",";
  print $sample;
}
print "\n";

foreach my $event(@splice_id_uniq)
{
  next unless scalar @{$count_spliced_event{$event}} >= 36;
  print $event;
  foreach my $sample (@samples_uniq){
    print ",";
    if($spliced_event{$event}{$sample} > 0 )
    {
      print $spliced_event{$event}{$sample};
    }
    else{
      print "0.001";
    }
  }
  print "\n";
  #print scalar @{$spliced_event{$event}},"*\n";
  #print join ",", @{$spliced_event{$event}},"\n";
}
__DATA__
print "SpliceID";
foreach my $sample(@sample_nonprog)
{
  print ",";
  print $sample;
}
foreach my $sample(@sample_prog)
{
  print ",";
  print $sample;
}
print "\n";
