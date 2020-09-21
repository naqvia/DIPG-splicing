#!/usr/bin/perl

##
## make_combined_tab.pl
##

my @ec_files= <*juncs.extracellular.wo.txt>;

foreach my $file (@ec_files)
{
  my $lib = $file;
  if($lib=~/(BS_[\w\d]+)/)
  {
    $lib = $1;
  }
  open(FIL,$file) || die("Cannot Open File");
  while(<FIL>)
  {
    chomp;
    print $lib,"\t",$_,"\n";
  }
}

__DATA__
PDX_1	chr1	159934424	159934425	ENSG00000085552.17:s:159934425-159934570_159934332-159934425	0.9103941917419434	-	chr1	159929761	159943149	1
