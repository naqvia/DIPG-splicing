#!/usr/bin/perl

#
# create_tab_for_splicing_classes.pl
#

my (%filter_mt,%filter_wt);

## add only dmg $samples
my $hist = "pbta-histologies.addedv16.dat";
open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
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

my (%file_paths_SE_controlBS,%file_paths_MXE_controlBS,%file_paths_A5SS_controlBS,%file_paths_A3SS_controlBS,%file_paths_RI_controlBS);
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
  $file_paths_SE_controlBS{$sample} = $file;

}
close(FIL);
open(FIL,"controlBS_vs_BS.MXE.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_MXE_controlBS{$sample} = $file;

}
close(FIL);
open(FIL,"controlBS_vs_BS.A5SS.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_A5SS_controlBS{$sample} = $file;

}
close(FIL);
open(FIL,"controlBS_vs_BS.A3SS.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_A3SS_controlBS{$sample} = $file;

}
close(FIL);
open(FIL,"controlBS_vs_BS.RI.rmat_files.txt") || die("Cannot Open File");
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
  $file_paths_RI_controlBS{$sample} = $file;

}
close(FIL);

my $A5SS_totals=0;
my $SE_totals=0;
my $MXE_totals=0;
my $RI_totals=0;
my $A3SS_totals=0;

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
  next unless ( $filter_mt{$sample} || $filter_wt{$sample});

  ## skip over G35 samples
  next if $filter_g35_hgg{$sample};

  ## control for HGGs non DMG are BR (brain homogenate)
  # if($filter_wt_hgg{$sample})
  # {
  #   $file = $file_paths_controlBR{$sample};
  #   #print "enter\t",$file,"\n";
  # }

  ## files for each type of splicing Type
  $file_SE  = $file_paths_SE_controlBS{$sample}; ## SE
  $file_MXE = $file_paths_MXE_controlBS{$sample};
  $file_RI  = $file_paths_RI_controlBS{$sample};
  $file_A5SS= $file_paths_A5SS_controlBS{$sample};
  $file_A3SS= $file_paths_A3SS_controlBS{$sample};

  my @files;
  push @files,$file_SE,$file_MXE,$file_RI,$file_A5SS,$file_A3SS;
  my $file_count = 1;

  foreach my $file(@files)
  {
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
        my $exon_up_start = $cols[3];
        my $exon_up_end    = $cols[4];
        my $exon_start    = $cols[5];
        my $exon_end      = $cols[6];
        my $exon_down_start = $cols[3];
        my $exon_down_end    = $cols[4];
        my $tumor_IJC     = $cols[14];
        my $tumor_SJC     = $cols[15];
        my $inc_from_len  = $cols[16];
        my $skip_from_len = $cols[17];


        my $pval     = $cols[-5];
        my $thr_diff = $cols[-1];

        ## only look at strong changes
        next unless (abs($thr_diff) >= .10);
        next unless ($pval<=0.05);

        ## remove last few items and print
        splice @cols, -5;


        ## 10 reads cut-off
        #next unless ($tumor_IJC >=10);
        #next unless ($tumor_SJC >=10);

        ##filter for microexons only
        #next if ($exon_end - $exon_start > 30);

        #print "*",$_,"\n";
        my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end."_".$exon_up_start."-".$exon_up_end."_".$exon_down_start."-".$exon_down_end."_".$str;
        $splice_id=~s/\"//g;


        if($file_count == 1)
        {
           $SE_totals++;
           #print $splice_id,"\tSE\n";
           splice @cols, -7;
           splice(@cols, 0, 1);
           print join "\t",@cols,"SingleE\n";
        }

        if($file_count == 2)
        {
           $MXE_totals++;
           #print $splice_id,"\tMXE\n";
           splice @cols, -7;
           splice(@cols, 0, 1);
           print join "\t",@cols,"MutualXE\n";
        }
        if($file_count == 3)
        {
           $RI_totals++;
           #print $splice_id,"\tRI\n";
           splice @cols, -7;
           splice(@cols, 0, 1);
           print join "\t",@cols,"ReIntr\n";
        }

        if($file_count == 4)
        {
           $A5SS_totals++;
           #print $splice_id,"\tA5SS\n";
           splice @cols, -7;
           splice(@cols, 0, 1);
           print join "\t",@cols,"A5SS\n";
        }

        if($file_count == 5)
        {
           $A3SS_totals++;
           #print $splice_id,"\tA3SS\n";
           splice @cols, -7;
           splice(@cols, 0, 1);
           print join "\t",@cols,"A3SS\n";
        }
      }
    }

    $file_count++;
  }
  close(SE_FIL);
}

#print "A5SS\tSE\tMXE\tA3SS\tRI\n";
#print $A5SS_totals,"\t",$SE_totals,"\t",$MXE_totals, "\t",$A3SS_totals, "\t",$RI_totals;
