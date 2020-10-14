#!/usr/bin/perl

#
# global_splicing_stats.pl
#

## tpm count table files
my ($polyA_tab, $stranded_tab) = ($ARGV[0], $ARGV[1]);

my @bs_sample_polya;
my @bs_sample_stranded;

my %sample_expr_totals;

## go through each tpm table and store tpm values for each gene
open(FIL, $polyA_tab) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  if($_=~/^BS/)
  {
    my @bs_ids = split;
    pop @bs_ids;
    push @bs_sample_polya,@bs_ids;
    #print "bs\t";
    foreach my $bs(@bs_ids)
    {
      #print "bs\t",$bs,"\n";
    }
  }

  my @cols = split;
  my $gene = pop @cols;
  #print "gene:".$gene,"\t";
  push @{$gene_tpms{$gene}},@cols;

  my $i = 0;
  foreach $tpm(@cols)
  {
    #print @bs_sample_polya[$i],":*",$tpm,",";
    if($tpm>1)
    {
        $sample_expr_totals{@bs_sample_polya[$i]}++;
    }
    $i++;
  }
  #print "\n";

}
close(FIL);

open(FIL, $stranded_tab) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  if($_=~/^BS/)
  {
    my @bs_ids = split;
    pop @bs_ids;
    push @bs_sample_stranded,@bs_ids;
    #print "bs\t";
    foreach my $bs(@bs_ids)
    {
      #print "bs\t",$bs,"\n";
    }
  }

  my @cols = split;
  my $gene = pop @cols;
  #print "gene:".$gene,"\t";
  push @{$gene_tpms{$gene}},@cols;

  my $i = 0;
  foreach $tpm(@cols)
  {
    #print @bs_sample_stranded[$i],":*",$tpm,",";
    if($tpm>1)
    {
        $sample_expr_totals{@bs_sample_stranded[$i]}++;
    }
    $i++;
  }
  #print "\n";


  push @{$gene_tpms{$gene}},@cols;
}
close(FIL);



## open histology file and annotate DIPG/DMG reads
my $hist = "pbta-histologies.addedv16.dat";
open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
  if($_=~/Diffuse\sintrinsic\spontine\sglioma/)
  {
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $ids_dmg{$id} = "H3K28";
    }
    else
    {
      $ids_dmg{$id} = "H3WT";
    }
  }
  else{
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $ids_dmg{$id} = "H3K28";
    }
    else
    {
      #$ids_dmg{$id} = "WT_H3";
    }
  }
}

## splicing files: psi and dPSI tsv files
#my @ctrl_bs_files     = <majiq_voila_tsv/dipg_deltapsinormal_control*tsv>;


my @psi_files_dmg         = <majiq_psi/dipg*psi.tsv>;
my @psi_files_hgg         = </Users/naqvia/Desktop/HGG/psi_majiq/*tsv>;
my @psi_files;

## store all PSI file names
push @psi_files,@psi_files_dmg;
push @psi_files,@psi_files_hgg;

#my @ctrl_bs_files     = <majiq_voila_tsv/dipg_deltapsinormal_control*tsv>;
#my @ctrl_bs_files = </Users/naqvia/Desktop/neoepitope_HGG/HGG-DIPG-neoepitope-non-denovo/control-BS/control-BS_vs_BS_*tsv>;
my @ctrl_bs_files = </Users/naqvia/Desktop/neoepitope_HGG/HGG-DIPG-neoepitope-non-denovo/control-BS/control-BS_vs_BS_*tsv>;

my %id_to_normal;
my %id_to_biospecimen;
my $hgg_dipg_manifest = $ARGV[2];

my @tsv_HGGs_DIPGs;
#push @tsv_HGGs_DIPGs,@ctrl_br_files;
push @tsv_HGGs_DIPGs,@ctrl_bs_files;
#push @tsv_HGGs_DIPGs,@ctrl_fb_files;
#push @tsv_HGGs_DIPGs,@ctrl_xx_files;
#push @tsv_HGGs_DIPGs,@srr4787052_files;
#push @tsv_HGGs_DIPGs,@nc_pub_files;




my %gene_aspliced;
my %asplicing_total;
my %splicing_total;
my %hgg_mappings;

open(MAP,"tsv_mappings.txt") || die("Cannot Open File");
while(<MAP>)
{
  chomp;
  ($name, $biospecimen) = split "\t";
  my $prefix = "\/Users\/naqvia\/Desktop\/HGG\/psi_majiq\/";
  my $full_name = $prefix."".$name;
  $full_name=~s/\s+//;
  $biospecimen=~s/\s+//;
  #print "full_name*".$full_name,"*\t*",$biospecimen,"*\n";
  $hgg_mappings{$full_name} = $biospecimen;
}

close(MAP);


## keep track and store all alternative spliced genes (not aberrantly but alternativly spliced genes)
foreach my $tsv(@psi_files)
{
	my @sep = split/\//,$tsv;
	my $biospecimen = "";

  if($tsv=~/(BS\_\w+)/)
    {
      $biospecimen = $1;
      $biospecimen=~s/BS_tumor_//;
      #print "bios:",$biospecimen,"\n";
    }
    else
    {
      #print "biosHGG:*",$tsv,"*\t*",$hgg_mappings{$tsv} ,"*\n";
      $biospecimen = $hgg_mappings{$tsv};
    }


  next unless $ids_dmg{$biospecimen};
  open(TSV,$tsv) || die("Cannot Open File");
  while(<TSV>)
  {
    chomp;
    next if ($_=~/\#/);

      my $gene = "";
      my @cols = split;
      if($_=~/(ENSG\d+)/)
      {
        #$gene = $cols[0];
        $gene = $1;
      }
      next if $gene_aspliced{$biospecimen}{$gene};
      $asplicing_total{$biospecimen}++;
      $gene_aspliced{$biospecimen}{$gene}=1;
    }
}

my %gene_spliced;
foreach my $tsv(@ctrl_bs_files)
{
	my @sep = split/\//,$tsv;
	my $biospecimen = $sep[-1];

  if($tsv=~/(BS\_\w+)/)
  {
    $biospecimen = $1;
    $biospecimen=~s/BS_vs_//;
  }
  else
  {
    #print "biosHGG:*",$tsv,"*\t*",$hgg_mappings{$tsv} ,"*\n";
    $biospecimen = $hgg_mappings{$tsv};
  }
	#my $normal_used_id = $id_to_normal{$biospecimen};
	#my $biospecimen    = $id_to_biospecimen{$biospecimen};

  #print $biospecimen,"*\n";
  next unless $ids_dmg{$biospecimen};
  #print "full_name*".$biospecimen,"*\n";
  print "*".$biospecimen,"\t",$tsv,"\n";
  open(TSV,$tsv) || die("Cannot Open File");
  while(<TSV>)
  {
    chomp;
    next if ($_=~/\#/);

      my $gene = "";
      my @cols = split;
      if($_=~/(ENSG\d+)/)
      {
        #$gene = $cols[0];
        $gene = $1;
      }
      #$gene = $cols[0];

      #print "sample ".$biospecimen." ".$gene,"\n";
      print "diff\t",$cols[0]."_".$cols[2],"\n";

      next if $gene_spliced{$biospecimen}{$gene};
      #print "sample ".$biospecimen." ".$gene,"\n";
      $splicing_total{$biospecimen}++;
      $gene_spliced{$biospecimen}{$gene}=1;
    }
}

print "Sample\tAlternative\tAberrant\n";
foreach my $sample(keys %sample_expr_totals)
{
  next unless $asplicing_total{$sample};
  next unless $ids_dmg{$sample};

  #print "totals\tsample:",$sample,"\texpr_totals:",$sample_expr_totals{$sample},"\taberrant:";
  #print $splicing_total{$sample},"\talternative:";
  #print $asplicing_total{$sample},"\t";

  print $sample,"\t";

  print sprintf ("%.2f",$asplicing_total{$sample}/$sample_expr_totals{$sample});
  print "\t";
  print sprintf ("%.2f",$splicing_total{$sample}/$asplicing_total{$sample});

  print "\n";
}

__DATA__
majiq_voila_tsv/dipg_deltapsinormal_control-BS_tumor_BS_A7Q8G0Y1.deltapsi.voila.tsv
