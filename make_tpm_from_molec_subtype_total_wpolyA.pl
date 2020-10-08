#!/usr/bin/perl
#use warnings;

my $tpm_1 = $ARGV[1];
my $tpm_2 = $ARGV[0];
my (@WT_samples, @MT_samples);
my %tpm;
my @total_samples;
my @total_samples_nonpoly;
my @total_samples_poly;

my (%mt_samples,%wt_samples);
my %filter_hgg_h3k8;
my @WT_samples_hgg;
my %filter_wt;

my (@nonpolyA_samples,@polyA_samples);

my @samples;
my %ids_dmg;

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
      push @genes, $gene;
      $i++;

    }
  }
}

##add only dmg $samples
my $hist = "pbta-histologies_w_added_col.tsv";
open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $id = $cols[0];

  my $polya_status = $cols[20];
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
      $status{$id} = "WT_hgg";
    }
  }
}

my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };
my @polyA_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @polyA_samples };
my @nonpolyA_samples_uniq = do { my %seen; grep { !$seen{$_}++ } @nonpolyA_samples };

my %poly_samples;
my %nonpoly_samples;

print "gene";
foreach my $polyA_sample(@polyA_samples_uniq)
{
  if($filter_wt_hgg{$polyA_sample})
  {
    #print ",polyA_hgg_wt:",$polyA_sample;
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
    #print ",nonpolyA_hgg_wt:",$nonpolyA_sample;
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







#die();

foreach my $gene(@genes_uniq)
{
    print $gene;

    foreach my $sample(@polyA_samples_uniq)
    {
      if($filter_wt_hgg{$sample})
      {
        #print ",";
        #my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        #print $rounded;
      }
    }

    foreach my $sample(@polyA_samples_uniq)
    {
      if($filter_wt{$sample})
      {
        print ",";
        my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        print $rounded;
      }
    }
    foreach my $sample(@polyA_samples_uniq)
    {
      if($filter_mt{$sample})
      {
        print ",";
        my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        print $rounded;
      }
    }

    ## non-polyA
    foreach my $sample(@nonpolyA_samples_uniq)
    {
      if($filter_wt_hgg{$sample})
      {
      #  print ",";
      #  my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      #  print $rounded; #int ",nonpolyA_hgg_wt:",$nonpolyA_sample;
      }
    }
    foreach my $sample(@nonpolyA_samples_uniq)
    {
      if($filter_wt{$sample})
      {
        print ",";
        my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        print $rounded;
      }
    }

    foreach my $sample(@nonpolyA_samples_uniq)
    {
      if($filter_mt{$sample})
      {
        print ",";
        my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        print $rounded;
      }
    }
    print "\n";

}

__DATA__
dipg_samples_w_molecsubtype.txt
