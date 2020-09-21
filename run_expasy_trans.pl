#!/usr/bin/perl

#
# run_expasy_trans.pl
#

my $file = $ARGV[0];

my $seqs = &get_fasta($file);
  foreach(keys %$seqs){
    print $_;
    print $$seqs{$_},"n";
  }



sub get_fasta{
  open(FILE, "<@_") or die("Cannot open FASTA file.n");
  my %seqs;
  my $header;
  my $first = 0;
  my @lines = <FILE>
  foreach my $line(@lines){
  chomp($line);
  if ($line =~ /^>/){
    $header = $line;
    $header =~ s/^>//;
    $header =~ s/s.*//;
    if ($first == 0){
      $first = 1;
    }
    next;
  }
  if ($first == 0){ die("Not a standard FASTA file.n"); }
    $seqs{$header} = $seqs{$header}.$line;
  }
  close(FILE);
  return %seqs;
}
