#!/usr/bin/perl

use POSIX;


while(<>)
{
  chomp;
  my @cols = split;
  my $gene = shift @cols;
  print $gene;
  foreach $expr(@cols)
  {
    print ",";
    my $ceil   = ceil($expr);
    print $ceil;
  }
  print "\n";
}

__DATA__
gene_counts_tpm.me_sfs.volplot.txt
