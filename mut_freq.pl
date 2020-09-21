#!/usr/bin/perl

open(FIL,"/Users/naqvia/Desktop/DIPG/mutation_freq_pedc.txt");
while(<FIL>)
{
  chomp;
  next if($_=~/Gene/);
  my ($gene,$freq,$cancer) = split;
  $freq=~s/\%//;

  my $splicing_freq = `grep \"^$gene\\s\" majiq_voila_tsv/dipg_deltapsinormal_control-BS_tumor_BS_\*tsv | awk '{print \$1}' | sort -u | wc -l`;
  $splicing_freq = sprintf ("%.2f",($splicing_freq/48)*100);

  if( ($splicing_freq >= 10) || ($freq >= 10)  )
  {
    print $gene,"\t";
    print $splicing_freq;
    print "\t";
    print $freq,"\t",$cancer,"\n";
  }
}
