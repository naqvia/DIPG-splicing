#!/usr/bin/perl

#
# report_DMG_splicing_events.pl
#

my $dir_file = $ARGV[0];

my %dmg_sample;

#add only DMG $samples
my $hist = "/Users/naqvia/Desktop/DIPG/pbta-histologies.addedv16.dat";

open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $id = $cols[0];
  if($_=~/Diffuse\sintrinsic\spontine\sglioma/)
  {
    my @cols = split "\t";
    my $id = $cols[0];


    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      #$filter_mt{$id} = "H3K28";
      #$status{$id} = "H3K28";
      $dmg_sample{$id} = $id;
    }
    else{
#      $filter_wt{$id} = "WT";
#      $status{$id} = "WT";
      $dmg_sample{$id} = $id;

    }
  }
  else{
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      #$filter_mt{$id} = "H3K28";
      #$status{$id} = "H3K28";
      $dmg_sample{$id} = $id;

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
  next unless $dmg_sample{$sample};
  print "sample\t";
  print $sample,"\n";
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
      next unless (abs($thr_diff) >= .05);
      next unless ($pval<=0.05);
      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);


      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end."_".$upstreamES."-".$upstreamEE."_".$downstreamES."-".$downstreamEE;
      $splice_id=~s/\"//g;
      print $splice_id,"\t";

      print $sample,"\t",$thr_diff,"\n";
      push @samples, $sample;
      push @splice_id, $splice_id;
      $samples_filter{$sample} = $sample;

      #print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";
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


__DATA__
