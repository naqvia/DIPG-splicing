#!/usr/bin/perl

#
# make_tab_for_total_diff_lsvs.pl
#

my $tsv_files = $ARGV[0];
my (%lsv_id_counts,%lsv_id_counts);

open(FIL,$tsv_files) || die("Cannot Open File");
while(<FIL>)
{
  my ($bs_id,$file) = split "\t";
  open(TSV,$file) || die("Cannot Open File");
  while(<TSV>)
  {
    chomp;
    next if($_=~/^LSV/);
    my @cols = split "\t";
    my $lsv_id = $cols[2];
    my $gene = $cols[0];
    my $lsv_id_v2 = $gene."\t".$lsv_id;
    $lsv_id_counts_v2{$lsv_id_v2}++;
    $lsv_id_counts{$lsv_id}++;
  }
  close(TSV);
}
close(FIL);

foreach my $lsv(keys %lsv_id_counts_v2)
{
  print $lsv,"\t",$lsv_id_counts_v2{$lsv},"\n";
  #print "total,".$lsv_id_counts{$lsv},"\n";

}
