#!/usr/bin/perl
#use warnings;

my $tpm_1 = $ARGV[1];
my $tpm_2 = $ARGV[0];
my (@WT_samples, @MT_samples);
my %tpm;
my @total_samples;
my (%mt_samples,%wt_samples);

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

print "gene";
foreach my $sample(@WT_samples)
{
  print ",";
  print $sample;
  $wt_samples{$sample} = 1;

}
foreach my $sample(@MT_samples)
{
  print ",";
  print $sample;
  $mt_samples{$sample} = 1;
}

foreach my $sample(@total_samples_uniq)
{
  next if $mt_samples{$sample};
  next if $wt_samples{$sample};
  print ",";
  print $sample,"";
}

print "\n";

foreach my $gene(@genes)
{
    print $gene;

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


    foreach my $sample(@total_samples_uniq)
    {
      next if $mt_samples{$sample};
      next if $wt_samples{$sample};
      print ",";
      my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      print $rounded;
    }

    print "\n";
  }

__DATA__
dipg_samples_w_molecsubtype.txt
