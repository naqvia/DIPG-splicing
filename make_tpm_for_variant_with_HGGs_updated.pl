#!/usr/bin/perl
#use warnings;

my $tpm_1 = $ARGV[1];
my $tpm_2 = $ARGV[0];
my (@WT_samples, @MT_samples);
my %tpm;
my @total_samples;
my (%mt_samples,%wt_samples);
my %filter_hgg_h3k8;
my @WT_samples_hgg;
my %filter_wt;


my @samples;
open(FIL1,$tpm_1) || die("Cannot Open File");
while(<FIL1>)
{
  chomp;

  if($_=~/^BS/)
  {
    @samples = split "\t";
    #print join "\t",@samples,"\n";
    pop @samples;
  }
  else
  {
    my @tpm = split "\t";
    my $gene = pop @tpm;
    my $i = 0;
    #print join "\t",@samples,"\n";
    foreach my $tpm (@tpm)
    {
      my $sample = $samples[$i];
      $tpm{$sample}{$gene} = $tpm;
      push @total_samples,$sample;
      #print $tpm,"\t",$sample,"\t",$gene,"\n";
      $i++;
    }
  }
}

my @samples;
open(FIL2,$tpm_2) || die("Cannot Open File");
while(<FIL2>)
{
  chomp;
  if($_=~/^BS/)
  {
    @samples = split "\t";
    pop @samples;
  }
  else
  {
    my @tpm = split "\t";
    my $gene = pop @tpm;
    my $i = 0;
    foreach my $tpm (@tpm)
    {
      my $sample = $samples[$i];
      $tpm{$sample}{$gene} = $tpm;
      push @total_samples,$sample;
      push @genes, $gene;
      $i++;

    }
  }
}

my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };
my @total_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @total_samples };



open(FILM,"dipg_samples_w_molecsubtype.txt") || die("Cannot Open File");
while (<FILM>)
{
  chomp;
  my ($bs_id,$pt_id,$status) = split "\t";
  if($status=~/K28/)
  {
    $status{$bs_id} = "H3K28M";
  #  print $bs_id," MT*\n";

    push @MT_samples, $bs_id;
  }
  else{
    $status{$bs_id} = "WT";
#    print $bs_id," WT*\n";

    push @WT_samples, $bs_id;
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
  #my ($HIST1H3B_ref,$HIST1H3B_alt =
  if( ($status_k28>=10))  #|| ($HIST1H3B<=10) || ($HIST1H3C<=10) || ($HIST2H3C<=10) || ($status_G35) <=10 || ($status_G35_2<=10 ))
  {
    #$status{$bs_id} = "H3K28M";
    #print $_,"\n";
    #push @MT_samples, $bs_id;
    #$filter_mt{$bs_id} = $bs_id;
    $filter_hgg_h3k8{$bs_id} = $bs_id;

  }
  else{
    $status{$bs_id} = "WT";
    push @WT_samples_hgg, $bs_id;
    $filter_wt{$bs_id} = $bs_id;

  }
}
close(FILM);


print "gene";
foreach my $sample(@WT_samples)
{
  print ",dmg_wt:";
  print $sample;
  $wt_samples{$sample} = 1;

}
foreach my $sample(@MT_samples)
{
  print ",dmg_h3k28:";
  print $sample;
  $mt_samples{$sample} = 1;
}

foreach my $sample(@total_samples_uniq)
{
  next if $mt_samples{$sample};
  next if $wt_samples{$sample};
  next if $filter_hgg_h3k8{$sample};
  print ",hgg_wt:";
  print $sample,"";
}

print "\n";
#die();

foreach my $gene(@genes)
{
    print $gene;

    foreach my $sample(@total_samples_uniq)
    {
      next if $mt_samples{$sample};
      next if $wt_samples{$sample};
      next if $filter_hgg_h3k8{$sample};

      print ",";
      my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      print $rounded;
    }

    foreach my $sample(@WT_samples)
    {
      print ",";
      #print $tpm{$sample}{$gene}

      my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      print $rounded;
    }

    foreach my $sample(@MT_samples)
    {
      print ",";
      my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      print $rounded;
    }




    print "\n";
  }

__DATA__
dipg_samples_w_molecsubtype.txt
