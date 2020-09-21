#!/usr/bin/perl

use warnings;
use strict;

while(<>)
{
	chomp;
	next if ($_=~/^#/);
	my @items = split "\t";
	my $deltapsi = $items[3];
	my @deltapsi_vals = split/\;/,$deltapsi;
	foreach my $deltapsi (@deltapsi_vals)
	{
		my $abs_deltapsi = abs($deltapsi);
		if($abs_deltapsi >= .30)
		{
			print $_,"\n";
		}
	}
}

__DATA__
cat dipg*deltapsi.voila.tsv
ITGA3	ENSG00000005884.17	ENSG00000005884.17:s:50086329-50087869	0.2516269147325176;-0.2516269227079554	0.9655640130491496;0.9655640401032901	1.2164119506101901e-05;1.2164845155426818e-05	0.694;0.306	0.977;0.023	s|1e1.2o2|1e2.1o1	False	False	True	2	3	1;1	chr17	+	50087869-50088225;50087869-50089110	50086329-50087869;50088213-50088366;50089110-50089815		http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A50089815-50086329