#!/usr/bin/perl

##
## make_psi_from_molec_subtype.pl
##

my $dir_file = $ARGV[0];

my (@MT_samples, @WT_samples);
my %tumor_IJCs;
my %tumor_SJCs;
my %splice_id_count;
my %filter_wt;
my %filter_mt;
my (%filter_hgg_h3k8,  %filter_g35_hgg, %filter_wt_hgg);
my %WT_samples_hgg;
my @g35_samples_hgg;

my (@polyA_samples, @nonpolyA_samples);
my (%file_paths_controlBR, %file_paths_controlBS);

##add only dmg $samples
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
my (%file_paths_controlBS,%file_paths_controlBR);

open(FIL,"controlBS_vs_BS.SE.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_controlBS{$sample} = $file;

}
close(FIL);

open(FIL,"controlBR_vs_BS.SE.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_controlBR{$sample} = $file;

}
close(FIL);

open(FIL,"controlBS_vs_BS.SE.rmat_files.txt") || die("Cannot Open File");
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

  ## skip unless H3K28, DMG, or HGG
  next unless ( $filter_wt_hgg{$sample} || $filter_mt{$sample} || $filter_wt{$sample});

  ## skip over G35 samples
  next if $filter_g35_hgg{$sample};

  ## control for HGGs non DMG are BR (brain homogenate)
  if($filter_wt_hgg{$sample})
  {
    $file = $file_paths_controlBR{$sample};
    #print "enter\t",$file,"\n";
  }

  #print "enter\t",$file,"\n";
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
      my $str  = $cols[4];

      my $exon_start    = $cols[5];
      my $exon_end      = $cols[6];
      my $tumor_IJC     = $cols[14];
      my $tumor_SJC     = $cols[15];
      my $inc_from_len  = $cols[16];
      my $skip_from_len = $cols[17];
      my $pval = $cols[18];
      my $thr_diff = $cols[-1];

      ## only look at strong changes
      next unless (abs($thr_diff) >= .10);
      next unless ($pval<=0.05);

      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);

      ##filter for microexons only
      next if ($exon_end - $exon_start > 30);

      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end;
      $splice_id=~s/\"//g;
      next unless ($pval<=0.05);
      push @samples, $sample;

      push @splice_id, $splice_id;
      $samples_filter{$sample} = $sample;

      #print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";

      push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      $splice_id_count{$splice_id}++;
      #$spliced_event{$splice_id}{$sample}
      #$tumor_IJCs{$splice_id}{$sample} = $thr_diff; #$thr_diff
      $tumor_IJCs{$splice_id}{$sample} = $thr_diff; #$thr_diff

      $tumor_SJCs{$splice_id}{$sample} = $tumor_SJC;
      $inc_from_lens{$splice_id} = $inc_from_len;
      $skip_from_lens{$splice_id} = $skip_from_len;
    }
  }
}
close(SE_FIL);

my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_id };
my @WT_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @WT_samples };
my @MT_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @MT_samples };
my @WT_samples_hgg_uniq = do { my %seen; grep { !$seen{$_}++ } @WT_samples_hgg };

## polyA vs non-polyA samples
my @polyA_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @polyA_samples };
my @nonpolyA_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @nonpolyA_samples };

## print headers for table
print "spliceID";
foreach my $polyA_sample(@polyA_samples_uniq)
{
  if($filter_wt_hgg{$polyA_sample})
  {
    print ",polyA_hgg_wt:",$polyA_sample;
  }
}
foreach my $polyA_sample(@polyA_samples_uniq)
{
  if($filter_wt{$polyA_sample})
  {
    print ",polyA_dmg_wt:",$polyA_sample;
  }
}
foreach my $polyA_sample(@polyA_samples_uniq)
{
  if($filter_mt{$polyA_sample})
  {
    print ",polyA_dmg_h3k28:",$polyA_sample;
  }
}

foreach my $nonpolyA_sample(@nonpolyA_samples_uniq)
{
  if($filter_wt_hgg{$nonpolyA_sample})
  {
    print ",nonpolyA_hgg_wt:",$nonpolyA_sample;
  }
}
foreach my $nonpolyA_sample(@nonpolyA_samples_uniq)
{
  if($filter_wt{$nonpolyA_sample})
  {
    print ",nonpolyA_dmg_wt:",$nonpolyA_sample;
  }
}

foreach my $nonpolyA_sample(@nonpolyA_samples_uniq)
{
  if($filter_mt{$nonpolyA_sample})
  {
    print ",nonpolyA_dmg_h3k28:",$nonpolyA_sample;
  }
}
print "\n";

foreach my $splice_id(@splice_id_uniq)
{
  next unless $splice_id_count{$splice_id} >=36;
  print $splice_id;

  foreach my $sample(@polyA_samples_uniq)
  {
    if($filter_wt_hgg{$sample})
    {
      print ",";
      if(abs($tumor_IJCs{$splice_id}{$sample} > 0) )
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0.001";
      }
    }
  }
  foreach my $sample(@polyA_samples_uniq)
  {
    if($filter_wt{$sample})
    {
      print ",";
      if(abs($tumor_IJCs{$splice_id}{$sample} > 0) )
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0.001";
      }
    }
  }
  foreach my $sample(@polyA_samples_uniq)
  {
    if($filter_mt{$sample})
    {
      print ",";
      if(abs($tumor_IJCs{$splice_id}{$sample} > 0) )
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0.001";
      }
    }
  }
  foreach my $sample(@nonpolyA_samples_uniq)
  {
    #print $sample
    if($filter_wt_hgg{$sample})
    {
      print ",";
      if(abs($tumor_IJCs{$splice_id}{$sample} > 0) )
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0.001";
      }
    }
  }
  foreach my $sample(@nonpolyA_samples_uniq)
  {
    if($filter_wt{$sample})
    {
      print ",";
      if(abs($tumor_IJCs{$splice_id}{$sample} > 0) )
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0.001";
      }
    }
  }
  foreach my $sample(@nonpolyA_samples_uniq)
  {
    if($filter_mt{$sample})
    {
      print ",";
      if(abs($tumor_IJCs{$splice_id}{$sample} > 0) )
      {
        print $tumor_IJCs{$splice_id}{$sample};
      }
      else{
        print "0.001";
      }
    }
  }
  print "\n";
}
