#!/usr/bin/perl

## process_for_SURVIV.pl

my $dir_file  = $ARGV[0];
my $histology = $ARGV[1];

my (@samples, @splice_id);
my %spliced_event;
my $samples_filter;
my (%sample_OS, %sample_OS_status);

my $tumor_IJCs;
my $tumor_SJCs;
my $inc_from_lens;
my $skip_from_lens;

my %filter_wt_hgg;

#add only HGG $samples
my $hist = "pbta-histologies.addedv16.dat";
open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $id = $cols[0];

  my $polya_status = $cols[21];
  if($polya_status=~/poly/)
  {
    push @polyA_samples, $id;
  }
  else
  {
    push @nonpolyA_samples, $id;
  }
  if($_=~/Diffuse\sintrinsic\spontine\sglioma/)
  {
    my @cols = split "\t";
    my $id = $cols[0];


    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $filter_mt{$id} = "H3K28";
      $status{$id} = "H3K28";
    }
    else{
      $filter_wt{$id} = "WT";
      $status{$id} = "WT";
    }
  }
  else{
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $filter_mt{$id} = "H3K28";
      $status{$id} = "H3K28";
    }
    elsif($molec_subtype=~/G35/)
    {
      $filter_g35_hgg{$id} = "G35";
      $status{$id} = "G35";
    }
    else
    {
      #print "*id:".$id,"\n";
      $filter_wt_hgg{$id} = "WT_H3";
      $status{$id} = "H3WT_HGG";
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
  next unless $filter_wt_hgg{$sample};

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
      my $upstreamES	  = $cols[7];
      my $upstreamEE	  = $cols[8];
      my $downstreamES	= $cols[9];
      my $downstreamEE  = $cols[10];


      my $tumor_IJC     = $cols[14];
      my $tumor_SJC     = $cols[15];
      my $inc_from_len  = $cols[16];
      my $skip_from_len = $cols[17];
      my $pval = $cols[18];
      my $thr_diff = $cols[-1];
      next unless (abs($thr_diff) >= .10);
      next unless ($pval<=0.05);
      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);


      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
      $splice_id=~s/\"//g;
      next unless ($pval<=0.05);
      push @samples, $sample;
      push @splice_id, $splice_id;
      $samples_filter{$sample} = $sample;

      print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";

      push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      $tumor_IJCs{$splice_id}{$sample} = $tumor_IJC;
      $tumor_SJCs{$splice_id}{$sample} = $tumor_SJC;
      $inc_from_lens{$splice_id} = $inc_from_len;
      $skip_from_lens{$splice_id} = $skip_from_len;
    }
  }
  close(SE_FIL);

  #print $file,"\n";
}
close(FIL);

open(FIL,$histology) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $biospecimen_id = $cols[0];
  my $OS = $cols[22];
  my $OS_status = $cols[23];
  next if ($OS=~/NA/);
 #print "hist *".$biospecimen_id,"*\t",$OS,"\t",$OS_status,"\n";

  if($samples_filter{$biospecimen_id})
  {
    $OS_days{$biospecimen_id}   = $OS;
    $OS_status{$biospecimen_id} = $OS_status;

    #print $biospecimen_id,"\t",$OS,"\t",$OS_status,"\n";
    $sample_OS{$biospecimen_id} = $OS;
    $sample_OS_status{$biospecimen_id} = $OS_status;
  }
}
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_id };

if($ARGV[2] == 1)
{
  print "PatientID\tTime\tEvent\n";

  foreach my $sample(@samples_uniq)
  {
    next unless $sample_OS{$sample}=~/\d+/;
    #next unless $OS_days{$sample};

    print $sample,"\t",$sample_OS{$sample},"\t";
    if($sample_OS_status{$sample} =~/LIVING/)
    {
      print "0\n";
    }
    else {
      print "1\n";
    }
  }
}
else{
  print "ID\tIJC\tSJC\tIncFormLen\tSkipFormLen\n";

  foreach my $splice_id (@splice_id_uniq)
  {
    print $splice_id, "\t";

    foreach my $sample(@samples_uniq)
    {
      next unless $OS_days{$sample}=~/\d+/;
      #print $sample,":";

      if($tumor_IJCs{$splice_id}{$sample} > 0)
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0";
      }
      print ",";

    }
    print "\t";

    foreach my $sample(@samples_uniq)
    {
      #print $sample,":";

      if($tumor_SJCs{$splice_id}{$sample} > 0)
      {
        print $tumor_SJCs{$splice_id}{$sample};
      }
      else{
        print "0";
      }
      print ",";
    }
    print "\t";

    print $inc_from_lens{$splice_id},"\t";
    print $skip_from_lens{$splice_id},"\n";
  }
}

__DATA__
For each sample-
    Open SE only events:
      —> Filter for cancer drivers
      —> Filter 10% diff. and strong P-vales (0.05)
      —> Each splicing event/gene
      —> Each patients against BR (significant splicing)

    Store:
    -Sample{splice_id} —> IJC
    -Sample{splice_id} —> SJC
    -Sample{splice_id} —> IncFormLen
    -Sample{splice_id} —> SkipFormLen

Open clinical file:
    -Sample —> survival status (dead/alive)
    -Sample —> OS
