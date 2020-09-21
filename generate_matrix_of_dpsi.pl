#!/usr/bin/perl

use strict;
use warnings;

#
# generate_matrix_of_dpsi.pl
#

my @tsv_files = <majiq_voila_tsv/*tsv>;

## hash data structure to store values
my %gene_lsv;
my %lsv_ctrl_psi;
my %lsv_case_psi;
my (@samples,@splicing_ids);
my %dpsi;

foreach my $file(@tsv_files)
{
  my $sample = "";
  if($file=~/(BS\_[\w+\d+]+)\.delta/)
  {
    #print "sample: ",$1,"\t";
    $sample = $1;
  }
  push @samples, $sample;

  open(FIL,$file) || die("Cannot Open File $file");
  while(<FIL>)
  {
      next if ($_=~/^#/);
      my (@content) = split "\t";
      my $gene_name = $content[0];
      my $lsv_coord = $content[2];
      my $psi_ctrl  = $content[6];
      my $psi_case  = $content[7];
      my $chr       = $content[15];
      my $strand    = $content[16];
      my ($ctrl_psi, $case_psi);
      my @psi_ctrl = split/;/,$psi_ctrl;
      my @psi_case = split/;/,$psi_case;

      #print $gene_name,"\t",$lsv_coord,"\t",$psi_ctrl,"\t",$psi_case,"\t";
      #print "\t";
      my $pos = 0;
      my $delta_psi = 0;
      while($pos < scalar(@psi_ctrl))
      {
        #print "enter\n";
        if($delta_psi > 0 || $delta_psi < 0)
        {
          if( abs($delta_psi) > abs($psi_ctrl[$pos] - $psi_case[$pos]) )
          {
            $delta_psi = $psi_ctrl[$pos] - $psi_case[$pos];
          }
        }
        else
        {
          $delta_psi = $psi_ctrl[$pos] - $psi_case[$pos];
        }
        $pos++;
      }
      #print "largest_dpsi: ".$delta_psi,"\t";
      my $first_event_dpsi = $psi_ctrl[0] - $psi_case[0];
      #print $first_event_dpsi,"\n";
      #print $gene_name,"_",$lsv_coord,"\t",$first_event_dpsi,"\n";
      my $splicing_id = $gene_name."_".$lsv_coord;
      push @splicing_ids, $splicing_id;
      $dpsi{$splicing_id}{$sample} = $first_event_dpsi;
  }
}

my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_ids_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_ids };

print "LSV\t";
foreach my $sample (@samples_uniq)
{
  print $sample,"\t";
}
print "\n";

#die();

foreach my $lsv (@splicing_ids_uniq)
{
  print $lsv,"\t";
  foreach my $sample (@samples_uniq)
  {
    if($dpsi{$lsv}{$sample})
    {
      print $dpsi{$lsv}{$sample};
    }
    else{ print "0";}
    #print $sample;
    print "\t";
  }
  print "\n";

}

__DATA__
#Gene Name	Gene ID	LSV ID	E(dPSI) per LSV junction	P(|dPSI|>=0.10) per LSV junction	P(|dPSI|<=0.05) per LSV junction	normal_control-BS E(PSI)	tumor_BS_B1C6GZ84 E(PSI)	LSV Type	A5SS	A3SS	ES	Num. Junctions	Num. Exons	De Novo Junctions	chr	strand	Junctions coords	Exons coords	IR coords	UCSC LSV Link
PKD1	ENSG00000008710.19	ENSG00000008710.19:t:2112788-2112963	0.8690560046520093;-0.8690815127683716	1.0;1.0	4.9333750481384e-05;4.933129482431742e-05	0.120;0.880	0.966;0.034	t|1e1.1o1|1e2.1o1	False	False	True	2	3	1;1	chr16	-	2112963-2113161;2112963-2113041	2112788-2112963;2113041-2113074;2113161-2113292		http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr16%3A2113292-2112788
