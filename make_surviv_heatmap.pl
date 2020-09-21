#!/usr/bin/perl

#
# make_surviv_heatmap.pl
#

my $samples = $ARGV[0];
my $output  = $ARGV[1];

my @sample;
my %pval_splicing;
my @splice_ids;

open(FIL,$samples) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my ($patient_id, $os, $status) = split "\t";
  push @samples, $patient_id;
}
close(FIL);

open(FIL,$output) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my (@cols) = split "\t";
  my $pval = $cols[-1];
  my $splice_id = $cols[0];
  my $calc_pval = log($pval) * -1;
  #print $_,"\n";
  #print "*",$calc_pval;
  push @{$pval_splicing{$splice_id}},$calc_pval;
  push @splice_ids, $splice_id;
}
close(FIL);
my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_ids };

print "SpliceID\t";
foreach my $sample(@samples)
{
  print $sample;
}

foreach my $event(@splice_id_uniq)
{
  print $event,"\t";
  print join "\t", @{$pval_splicing{$event}},"\n";
}
