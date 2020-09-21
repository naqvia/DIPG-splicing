#!/usr/bin/perl

## cross_check_in_DIPG_surv_exons.pl

my %sign_splice_ids;

my $dir_file = $ARGV[0];

open(INP, "/Users/naqvia/Desktop/SURVIV_analysis/SURVIV.out.BR.HGGs_nonDMG.SE.total.sign.txt") || die("Cannot Open File");
while(<INP>)
{
  chomp;
  my @cols = split "\t";
  my $splice_id = $cols[0];
  #print "SI*",$splice_id,"\n";
  $sign_splice_ids{$splice_id} = $splice_id;
}
open(FIL,$dir_file) || die("Cannot Open File $dir_file");
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
      my $exon_prev_start    = $cols[5];
      my $exon_prev_end      = $cols[6];
      my $exon_next_start    = $cols[5];
      my $exon_next_end      = $cols[6];

      my $tumor_IJC     = $cols[14];
      my $tumor_SJC     = $cols[15];
      my $inc_from_len  = $cols[16];
      my $skip_from_len = $cols[17];
      my $pval = $cols[18];
      my $thr_diff = $cols[-1];
      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end;
      $splice_id=~s/\"//g;

      next unless (abs($thr_diff) >= .10);
      next unless ($pval<=0.05);
      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);




      #next unless $sign_splice_ids{$splice_id};

      push @samples, $sample;
      push @splice_id, $splice_id;

      #print $_,"\n";
      print "*splice id: ",$splice_id,"\n";
       print $sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";

      push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      $tumor_IJCs{$splice_id}{$sample} = $tumor_IJC;
      $tumor_SJCs{$splice_id}{$sample} = $tumor_SJC;
      $inc_from_lens{$splice_id} = $inc_from_len;
      $skip_from_lens{$splice_id} = $skip_from_len;
    }
  }
}

  close(SE_FIL);
