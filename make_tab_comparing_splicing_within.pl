#!/usr/bin/perl

#
# make_tab_comparing_splicing_within.pl
#

my %filter_wt;
my %filter_mt;
my %filter_g35_hgg;
my %filter_wt_hgg;
my %splice_id_count;
my $h3k28_sample_count = 1;
my $h3wt_sample_count  = 1;
my %splice_id_count_mt;
my %splice_id_count_wt;

##add only dmg $samples
my $hist = "pbta-histologies.addedv16.dat";
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
      $filter_mt{$id} = "H3K28";
      $status{$id} = "H3K28";
      $h3k28_sample_count++;
    }
    else{
      $filter_wt{$id} = "WT";
      $status{$id} = "WT";
      $h3wt_sample_count++;
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
      $h3k28_sample_count++;

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
      $status{$id} = "WT_hgg";
    }
  }
}
my %file_paths_controlBS;
open(MAP,"controlBS_vs_BS.majiq_files.txt") || die("Cannot Open File");
while(<MAP>)
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
close(MAP);

open(FIL,"controlBS_vs_BS.majiq_files.txt") || die("Cannot Open File");
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

  ## skip over G35 samples
  next if $filter_g35_hgg{$sample};

  ## skip over HGGs H3WT
  next if $filter_wt_hgg{$sample};

  open(SE_FIL,$file);
  while(<SE_FIL>)
  {
    chomp;
    if($_=~/^\w/)
    {
      my @cols = split "\t";
      my $lsv = $cols[2];
      my $gene = $cols[0];
      my $lsv_id = $gene."_".$lsv;

      push @samples, $sample;
      push @splice_id, $lsv_id;

#     push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      if($filter_wt{$sample})
      {
        $splice_id_count_wt{$lsv_id}++;
        #print $sample,"\t",$splice_id,"\twt\n";
      }

      if($filter_mt{$sample})
      {
        $splice_id_count_mt{$lsv_id}++;
      }

      #$splice_id_count{$splice_id}{$sample}++;

    }
  }
}
close(SE_FIL);

my @samples_uniq   = do { my %seen; grep { !$seen{$_}++ } @samples   };
my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_id };


print "SpliceID,H3WT,H3K28\n";
foreach my $event(@splice_id_uniq)
{
  print $event;
  print ",";
  if($splice_id_count_wt{$event}) # %splice_id_count_mt
  {
    my $perc =  ($splice_id_count_wt{$event} / $h3wt_sample_count ) * 100;
    my $rounded = sprintf "%.0f", $perc;
    print $rounded;
    #print "(";
    #print $splice_id_count_wt{$event},",", $h3wt_sample_count;


  }
  else
  {
    print "0";
  }
  print ",";

  if($splice_id_count_mt{$event})
  {
    my $perc = ($splice_id_count_mt{$event} / $h3k28_sample_count ) * 100;
    my $rounded = sprintf "%.0f", $perc;
    print $rounded;
    #print ":";
    #print $splice_id_count_mt{$event},",", $h3k28_sample_count;
  }
  else
  {
    print "0";
  }
  print "\n";
}
