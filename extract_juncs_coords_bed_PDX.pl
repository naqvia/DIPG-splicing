#!/usr/bin/perl

#use warnings;
################################################################
## extract_juncs_coords_bed.pl
## 
## input: tsv MAJIQ file
###################################################################

#my @input = </Users/naqvia/Desktop/SJBALL/voila_non-denovo_non-ir_genome/neoepitopes/sjcloud_data/voila_p95/*deltapsi.txt>; #ProB_SJBALL020651.deltapsi_deltapsi.tsv
#my @input      = </Users/naqvia/Desktop/SJBALL/voila_non-denovo_non-ir_genome/neoepitopes/sjcloud_data/psi/voila_psi_di_dd_gencode/SJ*psi.txt>; #ProB_SJBALL020651.deltapsi_deltapsi.tsv

##ProB## 
my @input      = </Users/naqvia/Desktop/SJBALL/voila_non-denovo_non-ir_genome/neoepitopes/sjcloud_data/psi/voila_psi_di_dd_gencode/ProB*psi.txt>; #ProB_SJBALL020651.deltapsi_deltapsi.tsv
my @input = <*tsv>;

## go through each tsv file and parse
foreach my $file(@input)
{	
	my $lib = "";
	if($file=~/voila_psi_di_dd_gencode\/(SJ.+)\.psi/) ##Ctrl_SJBALL000013_D1.deltapsi.tsv
	{ 
		# print "enter ",$1,"\t";
		$lib = $1;
	}
	if($file=~/voila_psi_di_dd_gencode\/(ProB.+)\.psi/) ##Ctrl_SJBALL000013_D1.deltapsi.tsv
	{ 
		# print "enter ",$1,"\t";
		$lib = $1;
	}
	my $input_file = "/Users/naqvia/Desktop/SJBALL/voila_non-denovo_non-ir_genome/neoepitopes/sjcloud_data/psi/".$lib.".gencode_psi.thr01.bed";
	
	print "files\t",$input_file,"\t",$file,"\t",$lib,"\n";
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
		open(IN,">>",$input_file) || die("Cannot Open input file $input_file");
		my $i = 0;
		foreach my $junc(@junctions)
		{
			if($psi[$i]>=.10){
				my ($junc_start,$junc_end) = split/\-/,$junc;
				print IN $chr,"\t",$junc_start,"\t",$junc_end,"\t",$lsv_id,"\t",$psi[$i],"\t",$strand,"\n";
				$i++;
			}
			else{ $i++ }

		}
		#print IN $chr,"\t",$dpsi_exon_s,"\t",$dpsi_exon_e,"\t",$lsv_id,"\t",$dpsi[$index]+1,"\t",$strand,"\n";

		##die();
	}
}
__DATA__
		## to make sure exons are extracted
		if($strand=~/-/)
		{
			@lsv_type = reverse(@lsv_type);
		}


		foreach my $index(@dpsi_sign_indices)
		{
			#next unless($index=~/\d+/);
			#print "index_in_loop ",$index,"\t";

			#if($strand=~/-/){
			#	print "sign_section_-: ",$lsv_type[$index+1],"\t";
			#	my $section = $rev_lsv_type[$index-1];
			#}
			#else{
			my $section ="";
			if( ($lsv_id=~/\:t\:/) && ($strand=~/-/) ){
				$section = $lsv_type[$index];
				#print "section_t ",$section," from ",join "|",@lsv_type,"\t";
			}
			else{
				$section = $lsv_type[$index+1];
				#print "section ",$section," from ",$lsv_type,"\t";
			}
			my $exon_conn_index = "";
			if($section=~/\d+e(\d+)/)
			{

				$exon_conn_index = $1;
				#push @exon_conn,$exon_conn_index;
				## if target lsv then subtract one from exon connect because reference exon is the last exon
				if($lsv_type=~/t|/){ 
					$exon_conn_index= $exon_conn_index+1;
				}
				else{ 
					$exon_conn_index= $exon_conn_index;	
				}
				# print $lsv_id,"\tdpsis: ";
				# print join " ",@dpsi,"\t";
				# print "index_of_sign_dpsi: ",$index,"\t";
				# print "lsv_type: ",join "|",@lsv_type,"\t";
				# print "sign_section: ",$section,"\t";

				open(IN,">>",$input_file) || die("Cannot Open input file $input_file");
				if( ($strand=~/-/)&& ($lsv_id=~/\:s\:/) ){ 
					#print "exon_connect ".$exon_conn_index," ",$section,"\t",$strand,"\t",$exons[$exon_conn_index-1],"\n";
					my ($dpsi_exon_s, $dpsi_exon_e) = split/\-/,$exons[$exon_conn_index-1];
					next unless($dpsi_exon_s>=1);
					#print IN $chr,"\n";

					print IN $chr,"\t",$dpsi_exon_s,"\t",$dpsi_exon_e,"\t",$lsv_id,"\t",$dpsi[$index]+1,"\t",$strand,"\n";
					
				}
				elsif( ($strand=~/\+/) && ($lsv_id=~/\:s\:/) )
				{
					#print "exon_connect ".$exon_conn_index," ",$section,"\t",$strand,"\t",$exons[$exon_conn_index-1],"\n";
					my ($dpsi_exon_s, $dpsi_exon_e) = split/\-/,$exons[$exon_conn_index-1];
					next unless($dpsi_exon_s>=1);

					print IN $chr,"\t",$dpsi_exon_s,"\t",$dpsi_exon_e,"\t",$lsv_id,"\t",$dpsi[$index]+1,"\t",$strand,"\n";

				}
				elsif( ($strand=~/\-/) && ($lsv_id=~/\:t\:/) ){
					#print "exon_connect_neg ".$exon_conn_index," ",$section,"\t",$strand,"\t",$exons[$exon_conn_index-1],"\n";
					my ($dpsi_exon_s, $dpsi_exon_e) = split/\-/,$exons[$exon_conn_index-1];
					next unless($dpsi_exon_s>=1);
		
					print IN $chr,"\t",$dpsi_exon_s,"\t",$dpsi_exon_e,"\t",$lsv_id,"\t",$dpsi[$index]+1,"\t",$strand,"\n";
				}

				else{
					#print $lsv_id,"\t";
					#print "index_of_sign_dpsi: ",$index,"\t";
					#print "exon_connect ".$exon_conn_index," ",$section,"\t",$strand,"\t",$exons[$exon_conn_index-2],"\n";
					
					my ($dpsi_exon_s, $dpsi_exon_e) = split/\-/,$exons[$exon_conn_index-2];
					next unless($dpsi_exon_s>=1);

					print IN $chr,"\t",$dpsi_exon_s,"\t",$dpsi_exon_e,"\t",$lsv_id,"\t",$dpsi[$index]+1,"\t",$strand,"\n";


				}
				close(IN);
			}

		}

		#print join "\t",@exons," ec\n";

# 		foreach my $num(@exon_conn)
# 		{
			
# 			#print "exon_index ",$num, "\n"; 
# 			#print $exons[$num],"\n";
# 		}
# 		my @exon_conn = "";
# 		#print "\n";

# 		foreach $index(@dpsi_sign_indices)
# 		{
# 		#	print $lsv_id,"\t",$index,"\t";
# 		#	print $exons[$index],"\n";
# 		}
# 	#	print "\n";


# 		#foreach my $junc
# 		#print $gene,"\t",$ens,"\t",$lsv_id,"\n";

# 		# my ($Gene,$Gene_ID,$LSVID,$E,$P,$ProBE_1,$ProBE_2,$SJBALL013326,$LSVType,$A5SS,$A3SS,$ES,$NumJunctions,$NumExons,$DeNovoJunctions,$chr,$strand,$Junctionscoords,$Exonscoords,$IRcoords,$Voilalink) = split "\t";
# 		# my ($ens,$src,$coord) = split/\:/,$LSVID;
# 		# my $new_lsv = $Gene.":".$coord.":".$src;
# 		# my @exons = split/\;/,$Exonscoords;

# 		#print $Exonscoords,"*\n";
# 		my ($first_start, $first_end) = split/\-/, $exons[0];
# 		my ($last_start , $last_end ) = split/\-/, $exons[-1];

# 		my $lsv_start = $first_start;
# 		my $lsv_end   = $last_end;

# 		#print $lib,"\t",$chr,"\t",$lsv_start,"\t",$lsv_end,"\t",$lsv_id,"\t1\t",$strand,"\n";

# 		#print $chr,"\t",$lsv_start,"\t",$lsv_end,"\t", $Gene,"\t",$lsv_id,"\n";

# # 		foreach my $exon (@exons){
# # 			my ($start, $end) = split/\-/,$exon;
# # #			print $chr,"\t";
# # #			print $start,"\t",$end,"\t",$new_lsv,"\n";
# # 		}
	}
}
__DATA__
E(dPSI) per LSV junction tell you which path/junction to keep
ENSG00000012124:s:35820069-35820228
#Gene Name	Gene ID	LSV ID	E(dPSI) per LSV junction	P(|dPSI|>=0.20) per LSV junction	P(|dPSI|<=0.05) per LSV junction	Ctrl E(PSI)	SJBALL000013_D1 E(PSI)	LSV Type	A5SS	A3SS	ES	Num. Junctions	Num. Exons	De Novo Junctions	chr	strand	Junctions coords	Exons coords	IR coords	UCSC LSV Link
