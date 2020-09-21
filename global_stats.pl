#!/usr/bin/perl

#
# global_splicing_stats.pl
#

my ($polyA_tab, $stranded_tab) = ($ARGV[0], $ARGV[1]);
my @bs_sample_polya;
my @bs_sample_stranded;

my %sample_expr_totals;

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



my @ctrl_bs_files     = <majiq_voila_tsv/dipg_deltapsinormal_control*tsv>;
my @psi_files = <majiq_psi/dipg*psi.tsv>;

my %id_to_normal;
my %id_to_biospecimen;
my $hgg_dipg_manifest = $ARGV[2];

#push @tsv_HGGs_DIPGs,@ctrl_br_files;
push @tsv_HGGs_DIPGs,@ctrl_bs_files;
#push @tsv_HGGs_DIPGs,@ctrl_fb_files;
#push @tsv_HGGs_DIPGs,@ctrl_xx_files;
#push @tsv_HGGs_DIPGs,@srr4787052_files;
#push @tsv_HGGs_DIPGs,@nc_pub_files;

open(MAN,$hgg_dipg_manifest) || die("Cannot Open file $hgg_dipg_manifest");
while(<MAN>)
{
	chomp;
	next unless($_=~/tsv/);
  next unless ($_=~/control\-BS/);

	my @cols    = split/\,/;
	my $id          = $cols[0];
	my $name        = $cols[1];
	my $biospecimen = $cols[3];
	my $norm_id     = $cols[5];

  #print "name manifest ",$biospecimen,"\n";
	$id_to_normal{$biospecimen} = $norm_id;
	$id_to_biospecimen{$biospecimen} = $name;
	#$id_dipg{$biospecimen} = $biospecimen;
	push @biospecimen,$biospecimen;
}
close(MAN);


my %gene_aspliced;
my %asplicing_total;
my %splicing_total;

foreach my $tsv(@psi_files)
{
	my @sep = split/\//,$tsv;
	my $biospecimen = $sep[-1];
  if($tsv=~/(BS\_\w+)/)
  {
    $biospecimen = $1;
    $biospecimen=~s/BS_tumor_//;
  }
  #print "name ".$biospecimen,"\n";

	#my $normal_used_id = $id_to_normal{$name};
	#my $biospecimen    = $id_to_biospecimen{$name};

  #print "bio ",$biospecimen,"\n";
  open(TSV,$tsv) || die("Cannot Open File");
  while(<TSV>)
  {
    chomp;
    next if ($_=~/\#/);

      my $gene = "";
      my @cols = split;
      $gene = $cols[0];
      next if $gene_aspliced{$biospecimen}{$gene};
      #print "sample ".$biospecimen." ".$gene,"\n";
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
    $biospecimen=~s/BS_tumor_//;
  }
	#my $normal_used_id = $id_to_normal{$biospecimen};
	#my $biospecimen    = $id_to_biospecimen{$biospecimen};

  open(TSV,$tsv) || die("Cannot Open File");
  while(<TSV>)
  {
    chomp;
    next if ($_=~/\#/);

      my $gene = "";
      my @cols = split;
      $gene = $cols[0];
      #print "sample 1".$biospecimen." ".$gene,"\n";

      next if $gene_spliced{$biospecimen}{$gene};
      #print "sample ".$biospecimen." ".$gene,"\n";
      $splicing_total{$biospecimen}++;
      $gene_spliced{$biospecimen}{$gene}=1;
    }
}

foreach my $sample(keys %sample_expr_totals)
{
  next unless $asplicing_total{$sample};
  print "totals\t",$sample,"\t",$sample_expr_totals{$sample},"\t";
  print $splicing_total{$sample},"\t";
  print $asplicing_total{$sample},"\t";

  print sprintf ("%.2f",$asplicing_total{$sample}/$sample_expr_totals{$sample});
  print "\t";
  print sprintf ("%.2f",$splicing_total{$sample}/$asplicing_total{$sample});

  print "\n";
}

__DATA__
majiq_voila_tsv/dipg_deltapsinormal_control-BS_tumor_BS_A7Q8G0Y1.deltapsi.voila.tsv
