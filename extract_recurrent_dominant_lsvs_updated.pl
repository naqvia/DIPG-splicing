#!/usr/bin/perl

#use strict;
#use warnings;

#
# extract_recurrent_dominant_lsvs.pl
#

#my @tsv_files = <majiq_voila_tsv/*.tsv>;

## for denovo
#my @tsv_files = <majiq_denovo_voila_tsv/*.tsv>;

## for B-ALL
#my @tsv_files = <voila_deltapsi_disableIntron_disableDenovo_gcv31/ProB_*.tsv>;

#my @tsv_files = <voila_tsv/*non_denovo.lgg_wg_anno.voila_tsv.tsv>;

my $tsv_list = $ARGV[0];

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

my @tsv_files;
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
      $status{$id} = "WT_hgg";
    }
  }
}



my %file_paths_controlBS;
open(MAP,"controlBS_vs_BS.majiq_files.txt") || die("Cannot Open File");
while(<MAP>)
{
  chomp;
  my $file = $_;
  my $sample = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $sample = $1;
  }
  $file_paths_controlBS{$sample} = $file;

}
close(MAP);

my @ctrl_bs_files = </Users/naqvia/Desktop/neoepitope_HGG/HGG-DIPG-neoepitope-non-denovo/control-BS/control-BS_vs_BS_*tsv>;
foreach my $file(@ctrl_bs_files)
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
  next unless ($filter_wt{$biospecimen} || $filter_mt{$biospecimen});

  push @samples, $biospecimen;
  my $sample = $biospecimen;

  open(FIL,$file) || die("Cannot Open File $file");
  while(<FIL>)
  {
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

      print $gene_name,"\t",$lsv_coord,"\t",$psi_ctrl,"\t",$psi_case,"\t";
      print $exon_coords,"\t";
      #die;

      my $pos = 0;
      my $delta_psi = 0;
      my @relevant_psi_indices;
      my @relevant_dpsi;
      my %dpsi_index;

      ##look for PSI events that are 10% different
      while($pos < scalar(@psi_ctrl))
      {
        #print "enter\n";
        if( abs($psi_ctrl[$pos] - $psi_case[$pos]) >=.10 )
        {
          push @relevant_psi_indices, $pos;
          my $delta_psi = abs($psi_ctrl[$pos] - $psi_case[$pos]);
          push @relevant_dpsi, $delta_psi;
          $dpsi_index{$delta_psi} = $pos;
        }
        $pos++;
      }

      my @relevant_dpsi_sorted = sort { $a <=> $b } @relevant_dpsi;
      my $max_dpsi = $relevant_dpsi_sorted[-1];
      my $index_of_max_dpsi = $dpsi_index{$max_dpsi};

      print "max dpsi: ",$max_dpsi,"\tmax dpsi index: ",$index_of_max_dpsi,"\t";
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

  print "bed*".$chr{$event},"\t";
  if($event=~/(\d+\-\d+)\_(\d+\-\d+)/)
  {
    $lsv_coord = $1;
    $exon_coord = $2;
    #print "*".$exon_coord." ";
    my($start_exon_coord,$end_exon_coord) = split/\-/,$exon_coord;
    print $start_exon_coord,"\t",$end_exon_coord,"\t",$event,"\t",$avg_dpsi,"\t",$str{$event},"\n";
  }
}
## chr11	132335046	132335286	NTM_ENSG00000182667.14:s:132314552-132315070_132335046-132335286	0.640240740740741	+
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
