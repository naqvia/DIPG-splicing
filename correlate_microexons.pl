#!/usr/bin/perl

my $gene_expr = $ARGV[0];
my %tumor_expr;
my @genes;
my %normal_bs;
my %majiq_psi;
my %srrm4_expr;
my (%majiq_to_patient_id,%patient_to_majiq_id);
my %patient_rsem_id;
my %micro_genes;
my %ens_to_gene;
my %lsv_to_gene;
my @lsvs;

## map tpm ids and majiq id to patients
open(FIL, "/Users/naqvia/Desktop/DIPG/manifest/BSID_patientID_majiq_dpsi.txt") || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my ($majiq_id, $patient_id) = split;
  $majiq_to_patient_id{$majiq_id} = $patient_id;
  $patient_to_majiq_id{$patient_id} = $majiq_id;
  push @patient_ids, $patient_id;
}
close(FIL);

open(FIL, "/Users/naqvia/Desktop/DIPG/manifest/patientIDs_rsem.txt") || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my ($patient_id, $rsem_id) = split;
  $patient_rsem_id{$rsem_id} = $patient_id;
}
close(FIL);



open(FIL,$gene_expr) || die("Cannot Open File $gene_expr");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $normal_bs = $cols[6];
  my $gene = $cols[0];

  if($_=~/(ENSG\d+)\_/)

  {
    #$normal_bs{$gene} = $cols[6];
    push @genes, $gene;

    my $i = 0;
    my $ens = $1;

    foreach my $item(@cols)
    {
        my $name = $samples[$i];
        if($patient_rsem_id{$name})
        {
          print "i:".$i,"\n";
          my $rsem_name = $patient_rsem_id{$name};
          $tumor_expr{$ens}{$rsem_name} = $cols[$i];
          #print $gene,"\tsample name:",$name,"\t",$cols[$i],"\n";
        }
      $i++;
    }
  }
  elsif($_=~/^gene/){
    my @cols = split;
    shift @cols;
    #print join "\t",@cols,":cols\n";

    foreach my $sample(@cols)
    {
      #print $sample,"\n";

      push @samples,$sample;
    }
  }
  else{
    #print $_,"\n";
  }
}

## get psi values of BIN1
open(FIL,"/Users/naqvia/Desktop/DIPG/output/microexons_allpsi_patients.txt") || die("Cannot Open File");
while(<FIL>){
    chomp;
    my ($bs_id,$lsv,$psi) = split "\t";
    my (@events) = split/\;/,$psi;
    #print "psi_event:".$events[0],"\n";
    $majiq_psi{$bs_id}{$lsv} = $events[-1];
    push @lsvs,$lsv;
    if($lsv=~/(ENSG\d+)\./)
    {
      my $ens = $1;
      #print "ens".$ens,"\n";
      $lsv_to_gene{$lsv} = $ens;
      $micro_genes{$ens} = $ens;
    }

}

foreach my $lsv (@lsvs)
{
    my $out_file = "corr_v2/".$lsv.".corr_in.txt";
    open(OUT,">>".$out_file);

    #print "LSV ".$lsv,"\n";
    my $ens = "";
    if($lsv=~/(ENSG\d+)\./)
    {
      $ens = $1;
    }

    print $lsv,"\t",$ens,"*\n";
    print "ens".$ens,"\n";
        #print "ens".$ens,"\t",$i,"\n";
        my $i = 6;
        print $#samples," sample_counts\n";
        while($i<$#samples+1)
        {
          my $sample  = $samples[$i];
          my $patient_id =  $patient_rsem_id{$sample};
          my $majiq_id = $patient_to_majiq_id{$patient_id};
          my $psi = $majiq_psi{$majiq_id}{$lsv};
          my $expr = $srrm4_expr{$sample};
          #print "patient:",$patient_id,"\tmajiq:",$majiq_id,"\t",$expr,"\t",$psi,"\n";
          if($psi>.001){
            print OUT $lsv,"\t",$patient_id,"\t",$expr,"\t",$psi,"\n";
          }
      $i++;
      }
      close(OUT);
    }



__DATA__
foreach my $expr(@{$tumor_expr{"ENSG00000139767_SRRM4"}})
{
  my $sample  = $samples[$i];
  my $patient_id =  $patient_rsem_id{$sample};


  my $majiq_id = $patient_to_majiq_id{$patient_id};
  my $psi = $majiq_psi{$majiq_id};
  #next unless $majiq_psi{$majiq_id};

  if($psi>.001){
    print $patient_id,"\t",$expr,"\t",$psi,"\n";
  }
  $i++;
}

## store tpm values of SRRM4
foreach my $gene(@genes)
{
  my ($ens,$symb) = split/"_"/,$gene;

  if($gene=~/ENSG00000139767_SRRM4/)
  {
    #print $gene,"\n";

  #$ens_to_gene{$ens} = $gene;
  #if($micro_genes{$ens})
  #{
    my $i = 0;
    foreach my $expr(@{$tumor_expr{"ENSG00000139767"}})
    {
      my $sample = $samples[$i];
      #print $gene,"\t",$sample,"\t",$expr,"***\n";
      $srrm4_expr{$sample} = $expr;
      $i++;
    }
  }
}
