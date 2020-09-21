my %gene_counts;

while(<>)
{
  chomp;
  my @cols = split;
  my $lsv = shift @cols;
  my $norm_expr = shift @cols;
  next if($norm_expr == 0);
  print $lsv;
  foreach my $expr(@cols)
  {
    print "\t";
    #if($norm_expr == 0){ $norm_expr = 0.001; }
    print $expr/$norm_expr;
  }
  print "\n";
}

__DATA__
foreach my $gene(keys %gene_counts)
{
  print $gene,"\t",$gene_counts{$gene},"\n";
}

__DATA__
cat majiq_voila_tsv/*tsv | awk '{print $1"\t"$3}' | grep -v "#" | sort | uniq -c | sort -nr |

47 SMG7	ENSG00000116698.21:t:183548941-183549288
dpsi_SRSF.tab.v2.txt

SRSF10_ENSG00000188529.14:t:23974974-23975080	0	0	0	0	0	0	0	0.314	0	0	0	0	0	0	0	0	0	0	0	0	00	0	0	0	0	0	0.309	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.322	0	0	0
