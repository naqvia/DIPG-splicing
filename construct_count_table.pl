#!/usr/bin/perl

use strict;
use warnings;

my @tumor_files = <*genes.results>;
my @norm_files  = <normal_ctrl/*rsem.genes.results>;

my (@genes,@samples);
my (%gene_tpms, %ens,%gene);

foreach my $file (@norm_files)
{
    #chomp;
    #print $file,"**\t";
    my @file_path = split/\//,$file;
    my $sample = $file_path[-1];
    $sample =~s/\.rsem\.genes\.results//;
    $sample =~s/\.genes\.results//;
    push @samples,$sample;
    open(FIL,$file) || die("Cannot Open File $file");
    while(<FIL>){
        chomp;
        next if($_=~/^gene/);
        my ($gene,$transcript,$length,$effective_length,$expected_count,$TPM,$FPKM) = split "\t";
        if($gene=~/(ENSG\d+)/)
        {
            $gene = $1;
        }
        #print $sample,"\t",$gene,"\t",$TPM,"\n";
        push @genes, $gene;
        $gene_tpms{$sample}{$gene} = $TPM;
    }
}
foreach my $file (@tumor_files)
{
    #chomp;
    #print $file,"**\t";
    my @file_path = split/\//,$file;
    my $sample = $file_path[-1];
    $sample =~s/\.rsem\.genes\.results//;
    $sample =~s/\.genes\.results//;


    push @samples,$sample;
    open(FIL,$file) || die("Cannot Open File $file");
    while(<FIL>){
        chomp;
        next if($_=~/^gene/);
        my ($gene,$transcript,$length,$effective_length,$expected_count,$TPM,$FPKM) = split "\t";
        if($gene=~/(ENSG\d+)/)
        {
            $gene = $1;
        }
        #print $sample,"\t",$gene,"\t",$TPM,"\n";
        push @genes, $gene;
        $gene_tpms{$sample}{$gene} = $TPM;
    }
}

open(FIL,"ens_to_gene.dat");
while(<FIL>)
{
        chomp;
        my ($ens, $gene) = split;
        $ens{$ens} = $gene;
        $gene{$gene} = $ens;
}
close(FIL);


my @genes_uniq = do { my %seen; grep { !$seen{$_}++ } @genes };

print "gene\t";
print join "\t", @samples,"\n";
foreach my $gene (@genes_uniq)
{
  print $gene."_".$ens{$gene};
  foreach my $sample (@samples){
    print "\t";
    print $gene_tpms{$sample}{$gene};
  }
  print "\n";
}
