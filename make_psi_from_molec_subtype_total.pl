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
my (%filter_hgg_h3k8,  %filter_g35_hgg);
my %WT_samples_hgg;
my @g35_samples_hgg;

open(FILM,"dipg_samples_w_molecsubtype.txt") || die("Cannot Open File");
while (<FILM>)
{
  chomp;
  #print $_,"*\n";
  my ($bs_id,$pt_id,$status) = split "\t";
  if($status=~/K28/)
  {
    $status{$bs_id} = "H3K28M";
    #print $_,"\n";
    push @MT_samples, $bs_id;
    $filter_mt{$bs_id} = $bs_id;
  }
  else{
    $status{$bs_id} = "WT";
    push @WT_samples, $bs_id;
    $filter_wt{$bs_id} = $bs_id;

  }
}
close(FILM);

open(FILM,"input_dat/hgg_mpileup_hist_status.dat") || die("Cannot Open File");
while (<FILM>)
{
  chomp;
  next if $_=~/bs\-id/;

  #print $_,"*\n";
  my ($bs_id,$bam_path,$status_G35, $status_G35_2, $status_k28, $HIST1H3B, $HIST1H3C, $HIST2H3C) = split "\t";
  #my ($HIST1H3B_ref,$HIST1H3B_alt) = split "\/";

  if( ($status_k28>=10))  #|| ($HIST1H3B<=10) || ($HIST1H3C<=10) || ($HIST2H3C<=10) || ($status_G35) <=10 || ($status_G35_2<=10 ))
  {
    #$status{$bs_id} = "H3K28M";
    #print $_,"\n";
    #push @MT_samples, $bs_id;
    #$filter_mt{$bs_id} = $bs_id;
    $filter_hgg_h3k8{$bs_id} = $bs_id;

  }
  elsif($status_G35 >= 10 || $status_G35_2 >= 10)
  {
    $status{$bs_id} = "G35";
    push @G35_samples_hgg, $bs_id;
    $filter_g35_hgg{$bs_id} = $bs_id;
  }
  else{
    $status{$bs_id} = "WT";
    push @WT_samples_hgg, $bs_id;
    $filter_wt{$bs_id} = $bs_id;

  }
}
close(FILM);



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

  ##filter HGGs that have H3k28M mutation
  next if $filter_hgg_h3k8{$sample};

  #unless($filter_mt{$sample}){
  #  push @WT_samples,$sample;
  #}

  if ($filter_wt{$sample})
  {

  }
  elsif($filter_mt{$sample})
  {

  }
  elsif($filter_g35_hgg{$sample}){
    push @g35_samples_hgg,$sample
  }
  else{

    push @WT_samples_hgg,$sample;
  }
  #print "sample:",$sample,"\n";
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

      ## only look at strong changes
      next unless (abs($thr_diff) >= .10);
      next unless ($pval<=0.05);

      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);


      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end;
      $splice_id=~s/\"//g;
      next unless ($pval<=0.05);
      push @samples, $sample;

      push @splice_id, $splice_id;
      $samples_filter{$sample} = $sample;

    #  print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";

      push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      $splice_id_count{$splice_id}++;
      #$spliced_event{$splice_id}{$sample}
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

foreach my $sample(@WT_samples_uniq)
{
  print $sample,"*WT\n";
}
foreach my $sample(@MT_samples_uniq)
{
  print $sample,"*MT\n";
}
foreach my $sample(@WT_samples_hgg_uniq)
{
  print $sample,"*WT_hgg\n";
}
#die();

print "spliceID";
foreach my $sample(@WT_samples_hgg_uniq)
{
  print ",hgg_wt:",$sample;
}
foreach my $sample(@WT_samples_uniq)
{
  print ",dmg_wt:",$sample;
}
foreach my $sample(@MT_samples_uniq)
{
  print ",dmg_mt:",$sample;
}
foreach my $sample(@g35_samples_hgg)
{
  print ",hgg_g35:",$sample;
}
print "\n";


  foreach my $splice_id(@splice_id_uniq)
  {
      next unless $splice_id_count{$splice_id} >=36;
      print $splice_id;
      foreach my $sample(@WT_samples_hgg_uniq)
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

      foreach my $sample(@WT_samples_uniq)
      {
        #print $tpm{$sample}{$gene}

        print ",";
        if(abs($tumor_IJCs{$splice_id}{$sample}) > 0 )
        {
          print $tumor_IJCs{$splice_id}{$sample};
        }
        else{
          print "0.001";
        }

      }

      foreach my $sample(@MT_samples_uniq)
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

      foreach my $sample(@g35_samples_hgg)
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

      print "\n";

    }
