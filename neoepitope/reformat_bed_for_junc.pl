#!/usr/bin/perl

my @bed = <*.bed>;

foreach my $file(@bed)
{
	open (FIL, $file) || die("Cannot Open File $file");
	my $out_file = $file;

	$out_file =~s/bed/v2\.bed/; 
	print "out file".$out_file,"\n";
	open(IN,">>",$out_file) || die("Cannot Open input file $input_file");
	while(<FIL>)
	{
		chomp;
		my @cols = split "\t";
		my $chr  = $cols[0];
		my $start= $cols[1];
		my $end  = $cols[2];
		my $lsv  = $cols[3];
		my $psi  = $cols[4];
		my $str  = $cols[5];

		my $modified = $lsv; 
		$modified .= "_".$start."-".$end;
		if($str=~/\+/)
		{
			print IN $chr,"\t",$start,"\t",$start+1,"\t",$modified,"\t",$psi,"\t",$str,"\n";
		}
		else
		{
			print IN $chr,"\t",$end-1,"\t",$end,"\t",$modified,"\t",$psi,"\t",$str,"\n";	
		}
	}
	close($out_file);
}

__DATA__
chr17	63437560	63446245	ENSG00000008283.16:t:63437346-63437949	0.8174658417701721	-