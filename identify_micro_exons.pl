#!/usr/bin/perl

#
# identify_micro_exons.pl
#

my $input = $ARGV[0];

open(FIL, $input) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $exon_coords = $cols[5];

  my ($chr,$coord) = split/\:/,$exon_coords;
  my ($start_coord, $end_coord) = split/\-/,$coord;
  #print $chr,"\t",$start_coord,"\t",$end_coord,"\t";
  $exon_len = ($end_coord - $start_coord) + 1;
  next if($exon_len>30);
  if($exon_len % 3 > 0) {
    print "out_frame\t",$exon_len,"\t",$_,"\n";
  }
  else{
    print "in_frame\t",$exon_len,"\t",$_,"\n";
  }

}
__DATA__
ITGA3_ENSG00000005884.17:s:50086329-50087869_50089110-50089815   0.291354838709677       0.0091707810903135      3       31      chr17:50089110-50089815 +
