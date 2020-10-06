#!/usr/bin/perl

#use warnings;
##
## make_tab_for_uni_corr.pl
##

my $lsv_tab = $ARGV[0]; ## tab_global_splicing.total_lsv_names_updated.tsv
my $tsv_files = $ARGV[1]; ## dmg_file_names_total.dat

my %lsv_uni_filter;

open(FIL,$lsv_tab) || die("Cannot Open File $lsv_tab");
while(<FIL>)
{
  chomp;
  my ($gene, $lsv_id, $count) = split "\t";
  #my ($gene, $lsv_id, $count) = split(',', $_);
  if($count >= 53)
  {
    my $lsv_name = $gene."_".$lsv_id;
   #print "X".$count."\t".$lsv_name,"X\n";

    $lsv_uni_filter{$lsv_name} = $gene;
  }
}
close(FIL);

open(FIL,$tsv_files) || die("Cannot Open File");
while(<FIL>)
{
  my ($bs_id,$file) = split "\t";
  push @bs_ids,$bs_id;

  my $counter = 0;
  my %lsv_max_dpsi_index;

  open(TSV,$file) || die("Cannot Open File");
  while(<TSV>)
  {
    chomp;
    #next if($_=~/^LSV/);
    next if ($_=~/^#/);

    my @cols = split "\t";
    my $lsv_id = $cols[2];
    my $gene = $cols[0];
    my $lsv_name = $gene."_".$lsv_id;
    next unless $lsv_uni_filter{$lsv_name};
    push @total_lsv_names, $lsv_name;

    my (@content)  = split "\t";
    my $gene_name  = $content[0];
    my $lsv_coord  = $content[2];
    my $psi_ctrl   = $content[6];
    my $psi_case   = $content[7];
    my $num_exons  = $content[13];
    my $chr        = $content[15];
    my $strand     = $content[16];
    my $exon_coords= $content[18];

    my ($ctrl_psi, $case_psi);

    my @psi_ctrl = split/;/,$psi_ctrl;
    my @psi_case = split/;/,$psi_case;
    my @exon_coord = split/;/,$exon_coords;


    ## remove complex LSVs
    #next unless ( scalar(@psi_ctrl) == (scalar(@exon_coord)-1));

    #print $gene_name,"\t",$psi_ctrl,"\t",$psi_case,"\t";
    #print $exon_coords,"\t";

    my $pos = 0;
    my $delta_psi = 0;
    my @relevant_psi_indices;
    my @relevant_dpsi;
    my %dpsi_index;

    ##look for PSI events that are 10% different
    while($pos < scalar(@psi_ctrl))
    {

      push @relevant_psi_indices, $pos;
      my $delta_psi = abs($psi_ctrl[$pos] - $psi_case[$pos]);
      push @relevant_dpsi, $delta_psi;
      $dpsi_index{$delta_psi} = $pos;
      $pos++;
    }

    my @relevant_dpsi_sorted = sort { $a <=> $b } @relevant_dpsi;
    my $max_dpsi = $relevant_dpsi_sorted[-1];
    my $index_of_max_dpsi = $dpsi_index{$max_dpsi};

    #print $bs_id,"\t",$lsv_name,"\tmax dpsi: ",$max_dpsi,"\tmax dpsi index: ",$index_of_max_dpsi,"\tby_index: ";
    #print $relevant_dpsi[$index_of_max_dpsi],",\t";
    #print join "\t",@relevant_dpsi,"*\n";

  if($counter < 1)
  {
    $lsv_max_dpsi_index{$lsv_name} = $index_of_max_dpsi;
  }

    my $max_dpsi = $relevant_dpsi[$index_of_max_dpsi];
    $bs_splicing_psi{$bs_id}{$lsv_name} = $max_dpsi;
  }
}

my @total_lsv_names_uniq = do { my %seen; grep { !$seen{$_}++ } @total_lsv_names };

print "";
foreach my $lsv(@total_lsv_names_uniq)
{
  print ",";
  print $lsv;
}
print "\n";


foreach my $sample(@bs_ids)
{
  foreach my $lsv(@total_lsv_names_uniq)
  {

    print ",";
    print $bs_splicing_psi{$sample}{$lsv}
  }
  print "\n";
}


__DATA__
    #die;
    ## get relevant exon coordinates
    #     foreach my $psi_index(@relevant_psi_indices)
    #      {
      if( ($lsv_coord =~/\:s\:/) && ($strand=~/\+/) )
      {
        my $exon_coord_index = $index_of_max_dpsi + 1 ;
        print "enter max exon: ",$exon_coord[$exon_coord_index],"\n";
        my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
        $num_exons{$new_lsv_id} = $num_exons;
        $chr{$new_lsv_id} = $chr;
        $str{$new_lsv_id} = $strand;

        push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
               $exon_coord[$exon_coord_index];
        $splicing_event_count{$new_lsv_id}++;
        push @splicing_events, $new_lsv_id;
        push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;
      }

      # exons are in reverse order compared to psi vals
      elsif( ($lsv_coord =~/\:s\:/) && ($strand=~/\-/) ){
        my $exon_coord_index = $index_of_max_dpsi-1;
        my $num_exons = scalar(@exon_coord);
        my $num_psi   = scalar(@psi_ctrl);
        my $exon_coord_index = $num_exons - $num_psi;

        print "max exon: ",$exon_coord[$exon_coord_index],"\n";
        my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
        $num_exons{$new_lsv_id} = $num_exons;
        $chr{$new_lsv_id} = $chr;
        $str{$new_lsv_id} = $strand;

        push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
               $exon_coord[$exon_coord_index];
        $splicing_event_count{$new_lsv_id}++;
        push @splicing_events, $new_lsv_id;
        push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;

      }

      # exons are in reverse order compared to psi vals
      elsif( ($lsv_coord =~/\:t\:/) && ($strand=~/\-/) ){
        print scalar (@exon_coord)."&".scalar (@psi_ctrl),":";
        my $num_exons = scalar(@exon_coord);
        my $num_psi   = scalar(@psi_ctrl);
        my $exon_coord_index = $num_exons - $num_psi;
        print $exon_coord_index." ";
        #my $exon_coord_index = $index_of_max_dpsi;
        print "max exon: ",$exon_coord[$exon_coord_index],"\n";
        my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
        $num_exons{$new_lsv_id} = $num_exons;
        $chr{$new_lsv_id} = $chr;
        $str{$new_lsv_id} = $strand;

        push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
               $exon_coord[$exon_coord_index];
        $splicing_event_count{$new_lsv_id}++;
        push @splicing_events, $new_lsv_id;
        push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;

      }
      else{
        my $exon_coord_index = $index_of_max_dpsi;
        print "max exon: ",$exon_coord[$exon_coord_index],"\n";
        my $new_lsv_id = $gene_name."_".$lsv_coord."_".$exon_coord[$exon_coord_index];
        $num_exons{$new_lsv_id} = $num_exons;
        $chr{$new_lsv_id} = $chr;
        $str{$new_lsv_id} = $strand;

        push @{$relevant_splicing_events{$new_lsv_id}{$sample}},$max_dpsi,
               $exon_coord[$exon_coord_index];
        $splicing_event_count{$new_lsv_id}++;
        push @splicing_events, $new_lsv_id;
        push @{$splicing_event_deltapsis{$new_lsv_id}}, $max_dpsi;
      }
    #      }
    #die;
    print "largest_dpsi: ".$delta_psi,"\t";
  }
}
__DATA__

    #$bs_splicing_psi{$bs_id}{$lsv_name}
    $lsv_id_counts_v2{$lsv_id_v2}++;
    $lsv_id_counts{$lsv_id}++;
  }
  close(TSV);
}
close(FIL);

__DATA__
#Gene Name	Gene ID	LSV ID	E(dPSI) per LSV junction	P(|dPSI|>=0.10) per LSV junction	P(|dPSI|<=0.05) per LSV junction	control-BS_vs_ E(PSI)	BS_0ZA67BBC.540edb55-9185-4f8a-aaaf-fd59b465a15d.neoepitope_re-analysis.non_denovo E(PSI)	LSV Type	A5SS	A3SS	ES	Num. Junctions	Num. Exons	De Novo Junctions	chr	strand	Junctions coords	Exons coords	IR coords	UCSC LSV Link


ARHGEF12_120405933-120406141	PHLDB1_118639162-118639355	APMAP_24962925-24964022	DPF2_65343991-65344095	FLCN_17228026-17228161	AP1S2_15825806-15827372	ARHGEF10L_17632321-17632466	ARHGEF10L_17634835-17635016	ARVCF_19977415-19977586	ARVCF_19975686-19975762	TRIM9_50996064-50998188	TRIM9_50981800-50982311	NRXN2_64630402-64630573	NRXN2_64622753-64623996	CKAP5_46759269-46759442	SPTAN1_128591477-128591625	SPTAN1_128594122-128594373	ADD3_110133326-110135565	ITSN1_33799808-33799929	MADD_47308591-47308699	MADD_47324471-47325035	MADD_47309281-47309401	MADD_47326738-47327186	CALD1_134941092-134941237	REPS1_138921037-138921124	CLTB_176397607-176397718	CLTB_176392455-176392945	KDM1A_23059073-23059167	MEF2D_156475108-156475237	AP2A1_49802949-49803005	SBF1_50456492-50456673	CACNA1A_13212103-13212216	CACNA1A_13209312-13209498	ARHGEF7_111288354-111288443	CRTC1_18742910-18743026	CRTC1_18745823-18745960	SGCE_94585230-94585515	SHF_45175219-45175326	SHF_45172147-45172318	APBA2_29098490-29098576	CAMTA2_4969630-4969701	CAMTA2_4969150-4969415	SMG7_183545966-183546337	SMG7_183548941-183549288	CLTA_36211603-36212061	TMEM165_55423114-55424812	ABLIM1_114465465-114465827	ABLIM1_114451624-114451671	SUPT5H_39457675-39457949
0.99	0.236	0.998	0.103	0.994	0.838	0.018	0.999	0	0.986	0.036	0.996	0.497	0.505	0.998	0.091	0.917	0.949	0	0.156	0	0.888	0.996	0.959	0.618	0.074	0.885	0	0.933	0.823	0.081	0.09	0.855	0.073	0.082	0.985	0.965	0	0.818	0.973	0.012	0.936	0.251	0.755	0.972	0.04	0.325	0.809	0.98
