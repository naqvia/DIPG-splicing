#!/usr/bin/perl

my $gene_expr = $ARGV[0];

my (%majiq_to_patient_id,%patient_to_majiq_id);
my %patient_rsem_id;
my %micro_genes;
my %ens_to_gene;
my %lsv_to_gene;
my @lsvs;
my %rsem_to_patient_id;
my %tumor_expr;
my @genes;
my %ens_to_gene;
my  (@samples,@lsvs);
my (%H3K28_pt, %H3K28_bs);

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
  $rsem_to_patient_id{$patient_id} = $rsem_id;
}
close(FIL);


open(FIL,$gene_expr) || die("Cannot Open File $gene_expr");
while(<FIL>)
{
  chomp;
  my @cols = split;
  my $normal_bs = $cols[6];
  my $gene = $cols[0];

  push @genes, $gene;

  if($_=~/(ENSG\d+)\_/){
  my ($ensem,$official_name) = split/\_/,$gene;

    #$normal_bs{$gene} = $cols[6];

    my $i = 0;
    my $ens = $1;
    $ens_to_gene{$ens} = $official_name;

    foreach my $item(@cols)
    {
        my $name = $samples[$i];
        if($patient_rsem_id{$name})
        {
          my $rsem_name = $patient_rsem_id{$name};
          $tumor_expr{$rsem_name}{$gene} = $cols[$i];
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
      #print $sample,":sample\n";
      push @samples,$sample;
    }
  }
  else{
    #print $_,"\n";
  }
}

open(FIL,"dipg_samples_w_molecsubtype.txt");
while(<FIL>)
{
  chomp;
  my ($bs, $pt, $subtype) = split "\t";
  $H3K28_pt{$pt} = $bs;

  #print $bs,"\t",$pt."\t",$subtype,"***\n";
  if($subtype=~/H3\sK28/){
    $H3K28_bs{$bs} = $bs;
    #print $bs,"\t",$pt."\t",
    #$H3K28_bs{"H3K28"} = $bs;
  }
}



my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };

print "gene";
foreach my $sample(@samples_uniq)
{

  my $patient_id = $patient_rsem_id{$sample};
  #my $bs_sample  = $patient_to_majiq_id{$patient_id};
  my $bs_sample = $H3K28_pt{$patient_id};
  next unless ($bs_sample=~/\w+/);
  $bs_sample=~s/BS_tumor_//;
  unless ($H3K28_bs{$bs_sample})
  {

    print ",";
    print $bs_sample."_Non_H3K28";
  }
}
#print ",";

foreach my $sample(@samples_uniq)
{
  my $patient_id = $patient_rsem_id{$sample};
  my $bs_sample  = $patient_to_majiq_id{$patient_id};
  $bs_sample=~s/BS_tumor_//;
  next unless ($bs_sample=~/\w+/);

  if  ($H3K28_bs{$bs_sample})
  {
    print ",";
    print $bs_sample."_H3K28";
  }
}
print "\n";

foreach my $gene(@genes_uniq)
{
  next if($gene=~/gene/);

  print $gene;
  my %dupl;
  foreach my $sample(@samples_uniq)
  {
    next if($dupl{$bs_sample});
    my $patient_id = $patient_rsem_id{$sample};
    my $bs_sample  = $patient_to_majiq_id{$patient_id};
    next unless ($bs_sample=~/\w+/);

    unless  ($H3K28_bs{$bs_sample})
    {
      print ",";
      #print "gene:".$gene,"*\tsample:",$sample."\tpatient_id:".$patient_id,"*\tbs_id:".$bs_sample,"*\texpr:".$tumor_expr{$sample}{$gene},"\n";
      my $majiq_id = $patient_to_majiq_id{$sample};
      #print $patient,"\t";
        my $expr = $tumor_expr{$patient_id}{$gene};
        if($expr > 0)
        {
    #       print "gene:".$gene," ";
          print $expr;
        }
        else
        {
        #  print "sample:".$sample,"\t",$gene,"\n";
          print "0";
        }

    }
    $dupl{$bs_sample} = 1;

  }

  my %dupl = "";
  foreach my $sample(@samples_uniq)
  {
    my $patient_id = $patient_rsem_id{$sample};
    my $bs_sample  = $patient_to_majiq_id{$patient_id};
    next unless ($bs_sample=~/\w+/);
    next if($dupl{$bs_sample});

    if  ($H3K28_bs{$bs_sample})
    {
      print ",";
      #print "sample:",$sample.":".$lsv.":";  $tumor_expr{$rsem_name}{$gene}
      my $majiq_id = $patient_to_majiq_id{$sample};
      #print $patient,"\t";

        my $expr = $tumor_expr{$patient_id}{$gene};
        if($expr>0)
        {
    #       print "gene:".$gene," ";
          print $expr;
        }
        else
        {
          #print "gene:".$gene," ";
          print "0";
        }
    }
    $dupl{$bs_sample} = 1;
  }

  print "\n";
}

__DATA__




      if($psi_lsv{$lsv}{$sample} > 0){
        print $psi_lsv{$lsv}{$sample};
      }
      else{
        print "0";
      }
    }
  }
  foreach my $sample(@samples_uniq)
  {
    my $sample_to_check = $sample;
    $sample_to_check=~s/BS_tumor_//;
    if  ($H3K28_bs{$sample_to_check})
    {
      print ",";
      #print "sample:",$sample.":".$lsv.":";
      if($psi_lsv{$lsv}{$sample}){
        print $psi_lsv{$lsv}{$sample};
      }
      else{
        print "0";
      }
    }
  }
  print "\n";
}





foreach my $rsem_id (@samples_uniq)
{
    my $patient_id = $patient_rsem_id{$rsem_id};
    my $bs_id = $patient_to_majiq_id{$patient_id};

    unless($H3K28_bs{$bs})
    {
      print
    }

}

while(<>)
{
  chomp;
}

__DATA__
/Users/naqvia/Desktop/DIPG/gene_counts_tpm.txt
