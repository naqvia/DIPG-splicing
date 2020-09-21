#!/usr/bin/perl

use strict;
use warnings;

while(<>)
{
  chomp;
  if($_=~/LSV/) {
    print $_, "\n";
    next;
  }
  my @cols = split "\t";
  my $lsv = shift @cols;
  my $zero_counter = 0;
  foreach my $dpsi(@cols)
  {
    if($dpsi == 0){ $zero_counter++; }
  }

  if($zero_counter < 36)
  {
    print $_,"\n";
  }


}

__DATA__
EWSR1_ENSG00000182944.17:s:29273738-29273864	0.294	0.307	0.264	0.298	0.307	0.306	0.305	0.303	0	0.306	0.306	0.303	0.303	0.293	0.295	0.303	0.306	0.306	0.303	0.279	0.303	0.304	0.301	0.303	0.297	0.294	0.306	0.307	0.306	0.299	0.289	0.304	0.295	0.3	0.294	0.3	0.282	0.307	0.306	0.303	0.306	0.304	0.304	0.306	0.299	0.306	0.295
