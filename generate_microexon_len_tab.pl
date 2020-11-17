#!/usr/bin/perl

#
# generate_microexon_len_tab.pl
#

my $input_lsvs = $ARGV[0];
my %me_len_count;

my $me_outframe = 0;
my $me_inframe  = 0;

open(IN,$input_lsvs) || die("Cannot Open File");
while(<IN>)
{
  chomp;
  my @cols = split;
  my $me_len = $cols[-1];
  my $lsv_name = $cols[0];
  $me_len_count{$me_len}++;
  if( ($me_len % 3) == 0)
  {
    $me_inframe++;
  }
  else
  {
    $me_outframe++;
  }

}

print "Inframe\t";
print $me_inframe,"\n";

print "Outframe\t";
print $me_outframe,"\n";

my $len = 0;
while($len < 31)
{
  print $len,"\t";
  if($me_len_count{$len}){
    print $me_len_count{$len};
  }
  else{
    print "0";
  }
  print "\n";

  $len++;
}
__DATA__
PRKACB_ENSG00000142875.19:t:84179176-84179238_84175799-84175807	0.331875	0.00976699701900389	9	40	chr1:84175799-84175807	+	9
