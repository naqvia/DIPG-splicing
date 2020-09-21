#!/usr/bin/perl

##
## investigate_gene_expr.pl
##

my $gene_expr = $ARGV[0];
my %tumor_expr;
my @genes;
my %normal_bs;

open(FIL,$gene_expr) || die("Cannot Open File $gene_expr");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $normal_bs = $cols[6];
  my $gene = $cols[0];

  if($_=~/EN/)
  {
    $normal_bs{$gene} = $cols[6];
    push @genes, $gene;

    my $i = 0;
    foreach my $item(@cols)
    {
        if($i < 7)
        {
          $i++;
          next;
        }
        push @{$tumor_expr{$gene}}, $cols[$i];
        $i++;
    }
  }
}

foreach my $gene(@genes)
{
  print $gene,"\t",$normal_bs{$gene},"\t";
  my $avg_expr = &average(\@{$tumor_expr{$gene}});
  my $std_expr = &stdev(\@{$tumor_expr{$gene}});

  print $avg_expr,"\t",$std_expr,"\n";

  #print join "\t", @{$tumor_expr{$gene}},"\n";

  # my $sum = 0;
  # foreach my $expr(@{$tumor_expr{$gene}})
  # {
  #   $sum +=$expr;
  #
  # }
  # my $avg = $sum / scalar(@{$tumor_expr{$gene}});
}

## calculate stdev given an array of PSI values
sub stdev{
        my($data) = @_;
        if(@$data == 1){

                return "1";
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## calculate avg given an array of PSI values
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
