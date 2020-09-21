#!/usr/bin/perl

##
## make_final_tab_for_neojunc.pl
##

my $input = "DIPG_psi_neo_list_021720.txt";

## data structures
my @subtypes;
my %subtypes;
my (%gene_to_ens,%ens_to_gene);
my @ids;
my %entry_id;
my %total_psi;
my %total_num_junc;
my %tpms;
my %proB_tpms;

## get TPM values
my $file = "../gene_counts_tpm.txt";
open(FIL, $file) || die("Cannot Open TPM file");
while(<FIL>)
{
	chomp;
	if($_=~/ENS/){
		my @cols = split;
		my $gene_col = $cols[0];
		my ($ens,$gene) = split/\_/,$gene_col;


		$gene_to_ens{$gene} = $ens;
		$ens_to_gene{$ens} = $gene;
		my $i = 8;
		my $proB_tpms = "";
		foreach my $tpm(@cols)
		{
			#next if ($i == 0);
			#print $_,"\tTPM\t",$tpm,"\n";
			push @{$tpms{$ens}},$tpm;
			$i++;
			#last if($i == 39);
		}

		## store normal brain stem tpm
		push @{$proB_tpms{$ens}}, $cols[7];



		#$tpms_proB{$ens} = $proB_tpms;

	}


}



my %psi;

open(FIL,$input) || die("Cannot Open File $input");
while(<FIL>)
{
	chomp;
	my @cols = split "\t";
	my $lib        = $cols[0];
	my $chr        = $cols[1];
	my $start_junc = $cols[2];
	my $end_junc   = $cols[3];
	my $lsv        = $cols[4];
	my $psi        = $cols[5];
	my $str        = $cols[6];

	#next unless($uniq_lsvs{$lsv});
	my $id = $lsv."_".$chr.":".$start_junc."_".$end_junc;
	my $id = $lsv;
	my ($ens,$type,$coord) = split/\:/,$lsv;
	$ens =~s/\.\d+//;
	#print $gene,"\t",$lsv,"\t",$id,"\t",$
	push @ids, $id;

	$entry_id{$lib}{$id} = $_;
	$total_psi{$id}+=$psi;
	push @{$psi{$id}},$psi;
	$total_num_junc{$id}++;

}

my @uniq_ids    = do { my %seen; grep { !$seen{$_}++ } @ids };

foreach my $junc(@uniq_ids)
{
	#my $avg_psi = $total_psi{$junc} / $total_num_junc{$junc};

	my $avg_psi = &average(@psi{$junc});
	my $std_psi = &stdev(@psi{$junc});


	my (@splits) = split/\_/,$junc;
	my $lsv = $splits[0];
	my ($ens,$type,$coord) = split/\:/,$lsv;
	$ens =~s/\.\d+//;

	my $avg_tpm = &average(@tpms{$ens});
	my $std_tpm = &stdev(@tpms{$ens});
	my $avg_proB_tpms= "";
	my $avg_proB_tpms = &average(@proB_tpms{$ens});
	my $std_proB_tpms = &stdev(@proB_tpms{$ens});





	print $ens_to_gene{$ens},"\t",$lsv,"\t",$junc,"\t";
	print $avg_psi,"\t",$std_psi,"\t",$total_num_junc{$junc},"\t",$avg_tpm,"\t",$std_tpm,"\t",$avg_proB_tpms,"\t",$std_proB_tpms,"\n";
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                #die("Empty arrayn");
                $average = "NA*";
                return $average;
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

__DATA__
ENST00000374457	SRSF10	ENST00000374457	ENSG00000215699
SJBALL030295_D1	chr5	149782875	149784243	ENSG00000019582.14:s:149784243-149784330	0.8546985387802124	-	chr5	149781805	149786796	-	1368
