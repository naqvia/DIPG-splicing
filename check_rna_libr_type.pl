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
      push @total_samples_poly,$sample;
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
      push @total_samples_nonpoly,$sample;
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
  if($_=~/Diffuse\sintrinsic\spontine\sglioma/)
  {
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $ids_dmg{$id} = "H3K28";
    }
    else
    {
      $ids_dmg{$id} = "WT_H3";
    }
  }
  else{
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      #$ids_dmg{$id} = "H3K28";
    }
    else
    {
      #$ids_dmg{$id} = "WT_H3";
    }
  }
}

my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };
my @total_samples_poly_uniq = do { my %seen; grep { !$seen{$_}++ } @total_samples_poly };
my @total_samples_nonpoly_uniq = do { my %seen; grep { !$seen{$_}++ } @total_samples_nonpoly };

my %poly_samples;
my %nonpoly_samples;

print "gene";
foreach my $sample(@total_samples_poly_uniq)
{
  next unless $ids_dmg{$sample};

  print ",polyA:";

  print $sample;
  $poly_samples{$sample} = 1;

}
foreach my $sample(@total_samples_nonpoly_uniq)
{
  next unless $ids_dmg{$sample};

  print ",non_polyA:";

  print $sample;
  $nonpoly_samples{$sample} = 1;
}

print "\n";
#die();

foreach my $gene(@genes_uniq)
{
    print $gene;

    foreach my $sample(@total_samples_poly_uniq)
    {
      next unless $ids_dmg{$sample};

      #next if $nonpoly_samples{$sample};
      print ",";
      my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      print $rounded;
    }

    foreach my $sample(@total_samples_nonpoly_uniq)
    {
      next unless $ids_dmg{$sample};

      print ",";
      #print $tpm{$sample}{$gene}

      my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      print $rounded;
    }
    print "\n";
}

__DATA__
dipg_samples_w_molecsubtype.txt
