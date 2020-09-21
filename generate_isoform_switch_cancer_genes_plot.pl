#!/usr/bin/perl

#
# generate_isoform_switch_cancer_genes_plot.pl
#

my $isoform_switch_lsv = $ARGV[0];
my %dominant_lsvs;
my @tsv = <majiq_voila_tsv/dipg_deltapsinormal_control-BS_tumor_BS_*tsv>;
my @samples;
my (%dom_lsvs_in_sample,%lsvs_in_sample);

my $cancer_flag = $ARGV[1];
my (%ts_genes,%onco_genes);

  open(FIL1,"/Users/naqvia/Desktop/DIPG/onco_list.grep.txt")|| die("Cannot Open File");
  open(FIL2,"/Users/naqvia/Desktop/DIPG/ts_list.grep.txt")  || die("Cannot Open File");
  while(<FIL1>)
  {
    chomp;
    my $gene = $_;
    $gene=~s/^\^//;
    $gene=~s/\_$//;

    $onco_genes{$gene}=1;
  }
  close(FIL1);
  while(<FIL2>)
  {
    chomp;
    my $gene = $_;
    $gene=~s/^\^//;
    $gene=~s/\_$//;
    #print $gene,"\n";
    $ts_genes{$gene}=1;
  }
  close(FIL2);


open(FIL,$isoform_switch_lsv) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  $dominant_lsvs{$_} = 1;
}
close(FIL);

foreach my $tsv (@tsv)
{
  my $name = $tsv;
  push @samples,$name;
  open(FIL, $tsv)|| die("Cannot Open File");
  while(<FIL>)
  {
    chomp;
    my @cols = split "\t";
    my $lsv = $cols[2];
    my $gene = $cols[0];
    my $dpsi = $cols[3];
    my @dpsi_total = split/;/,$dpsi;

    my $thr_flag = 0;
    foreach my $dpsi(@dpsi_total)
    {
      if ( abs($dpsi) >=.30){
        $thr_flag = 1;
      }
    }
    next if($thr_flag == 0);
    #print $gene,"\n";
    if($ts_genes{$gene}) {
      $lsvs_in_sample{$name}++;
      if($dominant_lsvs{$lsv}){
        print "enter:",$gene,"\n";
        $dom_lsvs_in_sample{$name}++;
      }
    }

  }
close(FIL);
}
##sort by counts
foreach my $sample (sort { $lsvs_in_sample{$b} <=> $lsvs_in_sample{$a} } keys %lsvs_in_sample) {
    push @samples_ordered , $sample;
}



foreach my $sample(@samples_ordered)
{
  print $sample,"\t";
  my $non_dom_lsv_count = $lsvs_in_sample{$sample}-$dom_lsvs_in_sample{$sample};
  print $non_dom_lsv_count,"\tNon_switching\n";

  print $sample,"\t";
  print $dom_lsvs_in_sample{$sample},"\tSwitching\n";
}

__DATA__
ENSG00000001497.16:s:65524965-65525050
perl scripts/generate_isoform_switch_nums.pl dominant_events_lsvs.isoform_switch_lsv.txt
