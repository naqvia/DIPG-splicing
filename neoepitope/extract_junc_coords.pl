#!/usr/bin/perl

use strict;
use warnings;

##select which files you want to go through, non-denovo
my @tsv_input_deltapsi = </Users/naqvia/Desktop/DIPG/majiq_voila_tsv/dipg*voila.tsv>;
my %uniq_lsvs;

## go through each tsv file and parse
foreach my $file (@tsv_input_deltapsi)
{
	my $lib = $file;
	if($lib=~/tumor_(BS_[\w\d]+)/)
	{
		#print "yes";
		$lib = $1;
	}
	open(FIL,$file) || die("Cannot Open $file");
	while(<FIL>)
	{
		chomp;
		next unless($_=~/ENS/);
		my @cols = split "\t";
		#my $first = $cols[0];
		#my ($file,$gene) = split/\:/,$first;

		#$lib=$file;
		#print $file," file*\n";
		my $input_file = $lib.".dpsi_diffusage.bed";

		#my @cols = split "\t";
		#my $gene = $cols[0];BCAT1	ENSG00000060982.14	ENSG00000060982.14:t:24901814-24901885	0.4542948048741106;0.09948886318785066;-0.6446443923038088	0.9807393842981291;0.4923540397813762;0.9998449332649295	2.2594732947554125e-05;4.0262006062011246e-05;6.465083187690013e-06	0.167;0.095;0.738	0.650;0.251;0.102	t|1e1.1o1|1e2.1o1|1e3.1o1	False	False	True	3	4	1;1;1	chr12	-	24901885-24948927;24901885-24902947;24901885-24901963	24901814-24901885;24901963-24902388;24902947-24903075;24948927-24949459		http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr12%3A24949459-24901814
		my $ens = $cols[1];
		my $lsv_id=$cols[2];
		my $psi_per_junc = $cols[3];
		$psi_per_junc=~s/\s+//;

		my $ctrl_psi_per_junc = $cols[6];
		my $ball_psi_per_junc = $cols[7];
		$ctrl_psi_per_junc=~s/\s+//;
		$ball_psi_per_junc=~s/\s+//;

		my @ctrl_psi_per_junc = split/\;/,$ctrl_psi_per_junc;
		my @ball_psi_per_junc = split/\;/,$ball_psi_per_junc;

		#print $ctrl_psi_per_junc,"*\n";
		my @psi = split/\;/,$psi_per_junc;
		my @psi_sign_indices = ();
		my $i = 0;
		foreach my $junc(@psi)
		{

			push @psi_sign_indices,$i;
		 	$i++;
		}

		##identify exon coords that are affected
		#my (#Gene Name	Gene ID	LSV ID	E(dPSI) per LSV junction	P(|dPSI|>=0.10) per LSV junction	P(|dPSI|<=0.05) per LSV junction	ProB E(PSI)	PDX_29 E(PSI)	LSV Type	A5SS	A3SS	ES	Num. Junctions	Num. Exons	De Novo Junctions	chr	strand	Junctions coords	Exons coords	IR coords	UCSC LSV Link)
		my ($gene,$Gene_ID,$LSVID,$E,$V,$LSVType,$A5SS,$A3SS,$ES,$NumJunctions,$NumExons,$DeNovoJunctions,$chr,$strand,$Junctionscoords,$Exonscoords,$IRcoords) = split "\t";
		#my ($file,$gene) = split/\:/,$first_id;
		my @junctions = split/\;/,$Junctionscoords;
		my $strand = $cols[16];
		my $chr = $cols[15];

		##store junctions
		my @junctions = split/\;/,$cols[17];
		$ctrl_psi_per_junc=~s/\s+//;
		$ball_psi_per_junc=~s/\s+//;

		my @ctrl_psi_per_junc = split/\;/,$ctrl_psi_per_junc;


		my @lsv_type = split/\|/,$LSVType;


		print "lib ".$lib,"\tjuncs ";
		print join "\t",@junctions,"\t";
		print "\t",$lib,":",$input_file,"\n";

		#my $input_file = $file;
		#$input_file .= ".diff_usage.bed";

		open(IN,">>",$input_file) || die("Cannot Open input file $input_file");
		my $i = 0;
		foreach my $junc(@junctions)
		{

			if( ($ctrl_psi_per_junc[$i] < 0.01) && ($ball_psi_per_junc[$i] >= .10) )
			{
				my ($junc_start,$junc_end) = split/\-/,$junc;
				print IN $chr,"\t",$junc_start,"\t",$junc_end,"\t",$lsv_id,"_",$gene,"\t",$ball_psi_per_junc[$i]-$ctrl_psi_per_junc[$i],"\t",$strand,"\n";
				$i++;
			}

			elsif($ball_psi_per_junc[$i] >= $ctrl_psi_per_junc[$i]+.30)
			{
				my ($junc_start,$junc_end) = split/\-/,$junc;
				print IN $chr,"\t",$junc_start,"\t",$junc_end,"\t",$lsv_id,"_",,$gene,"\t",$ball_psi_per_junc[$i]-$ctrl_psi_per_junc[$i],"\t",$strand,"\n";
				$i++;
			}

			else{ $i++ }
		}
	}
}
__DATA__
## unique LSVs only in PDX -- compared to ProBs
my $uniq_to_PDX = "PDX_uniq_LSVs.txt";
open(FIL,$uniq_to_PDX) || die("Cannot Open $uniq_to_PDX");
while(<FIL>)
{
	chomp;
	my $lsv = $_;
	$lsv=~s/\s+//;
	#$uniq_lsvs{$lsv} = $lsv;
}

## go through each psi tsv file and parse
foreach my $file (@tsv_input_psi_BALL)
{
	my $lib = "";
	if($file=~/(PDX_.+)\.tsv/) ##Ctrl_SJBALL000013_D1.deltapsi.tsv
	{
		# print "enter ",$1,"\t";
		$lib = $1;
	}

	my $input_file = $lib.".psi_uniq_diffusage.bed";
	print $lib,"\n";

	open(FIL,$file) || die("Cannot Open $file");
	while(<FIL>)
	{
		chomp;
		#print $_,"\n";###################################################################
		next unless($_=~/ENS/);

		my @cols = split "\t";
		my $gene = $cols[0];
		my $ens = $cols[1];
		my $lsv_id=$cols[2];
		next unless $uniq_lsvs{$lsv_id};

		my $psi_per_junc = $cols[3];
		$psi_per_junc=~s/\s+//;
		#print $dpsi_per_junc,"*\n";
		my @psi = split/\;/,$psi_per_junc;
		my @psi_sign_indices = ();
		my $i = 0;
		foreach my $junc(@psi)
		 {

			push @psi_sign_indices,$i;
		 	$i++;
		 }

		##identify exon coords that are affected
		my ($Gene,$Gene_ID,$LSVID,$E,$V,$LSVType,$A5SS,$A3SS,$ES,$NumJunctions,$NumExons,$DeNovoJunctions,$chr,$strand,$Junctionscoords,$Exonscoords,$IRcoords) = split "\t";
		my @junctions = split/\;/,$Junctionscoords;
		my @lsv_type = split/\|/,$LSVType;

		print "lib ".$lib,"\tjuncs ";
		print join "\t",@junctions,"\t";
		print "\n";
		#open(IN,">>",$input_file) || die("Cannot Open input file $input_file");
		my $i = 0;
		foreach my $junc(@junctions)
		{
			if($psi[$i]>=.10){
				my ($junc_start,$junc_end) = split/\-/,$junc;
				print "IN chr",$chr,"\t",$junc_start,"\t",$junc_end,"\t",$lsv_id,"\t",$psi[$i],"\t",$strand,"\n";
				$i++;
			}
			else{ $i++ }

		}
		#print IN $chr,"\t",$dpsi_exon_s,"\t",$dpsi_exon_e,"\t",$lsv_id,"\t",$dpsi[$index]+1,"\t",$strand,"\n";

		##die();
	}
}


__DATA__



__DATA__

ProB E(PSI)	PDX_29 E(PSI)
0.918;0.082;0.000	0.578;0.420;0.002

	print "files\t",$input_file,"\t",$file,"\t",$lib,"\n";
	open(FIL,$file) || die("Cannot Open $file");
	while(<FIL>)
	{

		#Gene Name	Gene ID	LSV ID	E(dPSI) per LSV junction	P(|dPSI|>=0.10) per LSV junction	P(|dPSI|<=0.05) per LSV junction	ProB E(PSI)	PDX_29 E(PSI)	LSV Type	A5SS	A3SS	ES	Num. Junctions	Num. Exons	De Novo Junctions	chr	strand	Junctions coords	Exons coords	IR coords	UCSC LSV Link
