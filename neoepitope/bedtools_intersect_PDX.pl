#!/usr/bin/perl

## bedtools_intersect.pl

my @input_thr = <*.v2.bed>;
my $uniprot = "/Users/naqvia/Desktop/DIPG/unipLocExtra.hg38.col.txt";

foreach my $file(@input_thr)
{
	my $lib = $file;
	if($lib=~/(BS_[\w\d]+)/)
	{
		$lib = $1;
	}

	print "lib ".$lib,"*\n";
	#open(BED,$file) || die("Cannot Open bed file $bed");
	#while(<BED>)
	#{
		my $new_file = $lib.".dpsi.juncs.extracellular.wo.txt";
		#print $file,"\t",$new_file,"\n";
		print("bedtools intersect -wo -a $file -b $uniprot  > $new_file\n");
		system("bedtools intersect -wo -a $file -b $uniprot > $new_file");
	#}
}
