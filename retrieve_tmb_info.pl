#!/usr/bin/perl

#
# retrieve_tmb_info.pl
#

my $samples  = $ARGV[0];
my $tmb_file = $ARGV[1];

my @patient_ids;
my %patient_id_mappings;

open(FIL,$samples) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my ($bs, $patient_id) = split "\t";
  $patient_id_mappings{$patient_id} = $bs;
  push @patient_ids,$patient_id;
}
close(FIL);
open(FIL,$tmb_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my (@cols) = split "\t";
  my $bs = $cols[0];
  my $tmb = $cols[-1];
  $patient_id_mappings{$bs} = $tmb;
}
close(FIL);

my @patient_ids_uniq = do { my %seen; grep { !$seen{$_}++ } @patient_ids };

foreach my $patient_id(@patient_ids_uniq)
{

  my $bs = $patient_id_mappings{$patient_id};
  my $tmb = $patient_id_mappings{$bs};
  next unless $patient_id_mappings{$bs};
  print $patient_id,"\t";
  print $tmb,"\n";
}
