#!/usr/bin/perl

##
## identify_isoform_switch.pl
##

#!/usr/bin/perl

#use strict;
#use warnings;

#
# extract_recurrent_dominant_lsvs.pl
#

my @tsv_files = <majiq_voila_tsv/d*.tsv>; #BS_tumor_BS_0ZA67BBC.deltapsi


## hash data structure to store values
my %gene_lsv;
my %lsv_ctrl_psi;
my %lsv_case_psi;
my (@samples,@splicing_ids);
my %dpsi;

my (%splicing_event_count,%relevant_splicing_events);
my @splicing_events;
my %splicing_event_deltapsis;
my (%chr,%num_exons,%str);

my @tsv_files = </Users/naqvia/Desktop/neoepitope_HGG/HGG-DIPG-neoepitope-non-denovo/control-BS/control-BS_vs_BS_*tsv>;
my %ids_dmg;

## open histology file and annotate DIPG/DMG reads
##add only dmg $samples
my $hist = "pbta-histologies.addedv16.dat";
open(FIL,$hist) || die("Cannot Open File $hist");
while(<FIL>)
{
  chomp;
  my @cols = split "\t";
  my $id = $cols[0];

  my $polya_status = $cols[21];
  if($polya_status=~/poly/)
  {
    push @polyA_samples, $id;
  }
  else
  {
    push @nonpolyA_samples, $id;
  }
  if($_=~/Diffuse\sintrinsic\spontine\sglioma/)
  {
    my @cols = split "\t";
    my $id = $cols[0];


    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $filter_mt{$id} = "H3K28";
      $status{$id} = "H3K28";
    }
    else{
      $filter_wt{$id} = "WT";
      $status{$id} = "WT";
    }
  }
  else{
    my @cols = split "\t";
    my $id = $cols[0];
    my $molec_subtype = $cols[-1];
    if($molec_subtype=~/K28/)
    {
      $filter_mt{$id} = "H3K28";
      $status{$id} = "H3K28";
    }
    elsif($molec_subtype=~/G35/)
    {
      $filter_g35_hgg{$id} = "G35";
      $status{$id} = "G35";
    }
    else
    {
      #print "*id:".$id,"\n";
      $filter_wt_hgg{$id} = "WT_H3";
      $status{$id} = "H3WT_HGG";
    }
  }
}

foreach my $file(@tsv_files)
{

  my @sep = split/\//,$file;
  my $biospecimen = $sep[-1];


  if($file=~/(BS\_\w+)/)
  {
    $biospecimen = $1;
    $biospecimen=~s/BS_vs_//;
  }
  else
  {
    #print "biosHGG:*",$tsv,"*\t*",$hgg_mappings{$tsv} ,"*\n";
    $biospecimen = $hgg_mappings{$file};
  }

  #print "biospec:".$biospecimen,"\n";

  next unless ($filter_mt{$biospecimen} || $filter_wt{$biospecimen});

  push @samples, $biospecimen;
  my $sample = $biospecimen;

  open(FIL,$file) || die("Cannot Open File $file");
  while(<FIL>)
  {
    chomp;
      next if ($_=~/^#/);
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
      next unless ( scalar(@psi_ctrl) == (scalar(@exon_coord)-1));

      my $pos = 0;
      my $delta_psi = 0;
      my @relevant_psi_indices;
      my @relevant_dpsi;
      my %dpsi_index;

      ## find biggest PSI in normal and tumor @samples
      my @psi_sorted_ctrl = sort { $a <=> $b } @psi_ctrl;
      my $max_psi_ctrl    = @psi_sorted_ctrl[-1];

      my @psi_sorted_tumor = sort { $a <=> $b } @psi_case;
      my $max_psi_tumor    = @psi_sorted_tumor[-1];

      my $i = 0;
      foreach my $psi(@psi_ctrl)
      {
          if($psi=~/$max_psi_ctrl/)
          {
            $max_psi_ctrl_index = $i;
          }
          else
          {
            $i++;
          }
      }

      my $j = 0;
      foreach my $psi(@psi_case)
      {
          if($psi=~/$max_psi_tumor/)
          {
            $max_psi_tumor_index = $j;
          }
          else
          {
            $j++;
          }
      }


      #go to next event if the same exons are dominant in both conditions
      if ($max_psi_tumor_index == $max_psi_ctrl_index){
        #print $sample,"\t",$lsv_coord,"\tSame\n";
        #print $_,"\n";
      }
      else{
        #print $sample,"\t",$lsv_coord,"\tSwitch\n";
        print $_,"\n";
        #print $gene_name,"\t",$lsv_coord,"\t",$psi_ctrl,"\t",$psi_case,"\t";
        #print $exon_coords,"\t";
        #print $max_psi_ctrl_index,",",$max_psi_tumor_index,"\t";
        #print "*",$lsv_coord,"\n";
      #die;
    }

  }
}

__DATA__
perl scripts/identify_isoform_switch.pl  |  grep "*" | perl -pe 's/\*//'  | sort -u > dominant_events_lsvs.isoform_switch_lsv.txt


## make arrays unique
my @samples_uniq = do { my %seen; grep { !$seen{$_}++ } @samples };
my @splicing_events_uniq = do { my %seen; grep { !$seen{$_}++ } @splicing_events };

foreach my $event (@splicing_events_uniq)
{
  ##compute avg and stdevs of each lsv dPSI
  print "event: ",$event,"\t";
  print join "\t",@{$splicing_event_deltapsis{$event}},"\n";
  my $avg_dpsi = &average(\@{$splicing_event_deltapsis{$event}});
	my $std_dpsi = &stdev(\@{$splicing_event_deltapsis{$event}});

  #next unless $splicing_event_count{$event} < 9;
  print "*",$event,"\t",$avg_dpsi,"\t",$std_dpsi,"\t",$num_exons{$event},"\t";
  print $splicing_event_count{$event};
  my ($gene,$lsv,$exon_coord) = split/\_/,$event;
  print "\t",$chr{$event},":",$exon_coord,"\t",$str{$event},"\n";

}

## calculate stdev given an array of PSI values
sub stdev{
        my($data) = @_;
        if(@$data == 1){

                return "1";
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## calculate avg given an array of PSI values
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
__DATA__
perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | awk -F "\t" '{if ($5 == 47){ print $0 }}' | awk '{print $6"\t"$1"\t"$2"\t"$7}'| perl -pe 's/(chr[\d+XY]+)\:/$1\t/g' | awk -F "-" '{print $1"\t"$2"-"$3"-"$4}' | awk '{if($6 == "+"){ print $0 } else {print $0"-" }}' | perl -pe 's/\*//' >
