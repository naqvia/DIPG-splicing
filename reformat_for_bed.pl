#!/usr/bin/perl

while(<>)
{
  chomp;
  my @items = split "\t";
  my $coord_id = $items[5];
  my $id     = $items[0];
  my $avg_dpsi = $items[1];
  my $str = $items[6];

  ## chr17:50089110-50089815
  my ($chr,$coord) = split/\:/,$coord_id;
  my ($start_coord, $end_coord) = split/\-/,$coord;
  print $chr,"\t",$start_coord,"\t",$end_coord,"\t";

  $id =~s/\*//;
  print $id,"\t", $avg_dpsi,"\t",$str,"\n";

}

__DATA__
*ITGA3_ENSG00000005884.17:s:50086329-50087869_50089110-50089815  0.291354838709677       0.0091707810903135      3       31      chr17:50089110-50089815 +
