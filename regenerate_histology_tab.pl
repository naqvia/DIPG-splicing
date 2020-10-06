#!/usr/bin/perl

##
## regenerate_histology_tab.pl
##

my $pbta_hist = $ARGV[0];

my (%dmg_samples, %hgg_samples);
my %total_samples;
my (%filter_hgg_wt, %filter_hgg_g35, %filter_hgg_g35, %filter_mt, %filter_wt);
open(FILM,"dipg_samples_w_molecsubtype.txt") || die("Cannot Open File");
while (<FILM>)
{
  chomp;
  #print $_,"*\n";
  my ($bs_id,$pt_id,$status) = split "\t";
  $total_samples{$bs_id} = $bs_id;
  if($status=~/K28/)
  {
    $status{$bs_id} = "H3K28M";
    #print $_,"\n";
    push @MT_samples, $bs_id;
    $filter_mt{$bs_id} = $bs_id;
  }
  else{
    $status{$bs_id} = "WT";
    push @WT_samples, $bs_id;
    $filter_wt{$bs_id} = $bs_id;

  }
}
close(FILM);

open(FILM,"input_dat/hgg_mpileup_hist_status.dat") || die("Cannot Open File");
while (<FILM>)
{
  chomp;
  next if $_=~/bs\-id/;

  #print $_,"*\n";
  my ($bs_id,$bam_path,$status_G35, $status_G35_2, $status_k28, $HIST1H3B, $HIST1H3C, $HIST2H3C) = split "\t";
  #my ($HIST1H3B_ref,$HIST1H3B_alt) = split "\/";
  $total_samples{$bs_id} = $bs_id;

  if( ($status_k28>=10))  #|| ($HIST1H3B<=10) || ($HIST1H3C<=10) || ($HIST2H3C<=10) || ($status_G35) <=10 || ($status_G35_2<=10 ))
  {
    $filter_hgg_h3k8{$bs_id} = $bs_id;
  }
  elsif($status_G35 >= 10 || $status_G35_2 >= 10)
  {
    $filter_hgg_g35{$bs_id} = $bs_id;
  }
  else{
    $filter_wt_hgg{$bs_id} = $bs_id;
  }
}
close(FILM);

open(WGS_HGG,"hgg_non-dmg_samples.RNAseq_and_patientIDs_withWGS.txt") || die("Cannot Open File");
while (<WGS_HGG>)
{
  chomp;
  my ($bs_id,$patient_id) = split "\t";
  $total_samples{$bs_id} = $bs_id;
  $filter_wt_hgg{$bs_id} = $bs_id;
}

open(HIST,$pbta_hist) || die("Cannot Open File $pbta_hist");
while(<HIST>)
{
  chomp;
  my @cols = split "\t";
  my $kids_first = $cols[0];
  my $bs_id      = $cols[1];
  next unless $_=~/RNA\-Seq/;
  #next unless $_=~/$kids_first/;
  next unless $total_samples{$kids_first};

  if($filter_wt{$kids_first})
  {
    print $_,"\t*H3-WT\n";
  }
  elsif($filter_mt{$kids_first})
  {
    print $_,"\t*H3K28\n";
  }
  elsif($filter_hgg_h3k8{$kids_first})
  {
    print $_,"\t*H3K28\n";

  }
  elsif($filter_hgg_g35{$kids_first})
  {
    print $_,"\t*G35\n";

  }
  elsif($filter_wt_hgg{$kids_first})
  {
    print $_,"\t*H3-WT\n";

  }
  else{
    #print "*",$_,"\n";
  }
}


__DATA__
Kids_First_Biospecimen_ID	sample_id	aliquot_id	Kids_First_Participant_ID	experimental_strategy	sample_type	composition	tumor_descriptor	primary_site	reported_gender	race	ethnicity	age_at_diagnosis_days	disease_type_old	disease_type_new	short_histology	broad_histology	broad_composition	Notes	germline_sex_estimate	RNA_library	OS_days	OS_status	cohort	age_last_update_days	source_text_tumor_descriptor	cancer_predispositions	seq_center	normal_fraction	tumor_fraction	glioma_brain_region	tumor_ploidy	molecular_subtype
