#!/usr/bin/perl

## make_tab_for_pegasus.pl

my ($samples_file, $expr_file,$pathway_file) = ($ARGV[0], $ARGV[1], $ARGV[2]);
my @pre_dipg_samples;
my (@dipg_samples,@hgg_samples,@lgg_samples);

## files from MAJIQ output for each cancer Type

my @tsv_majiq_dipg = </Users/naqvia/Desktop/DIPG/majiq_voila_tsv/dipg_deltapsinormal_control-BS_tumor_BS_0ZA67BBC.deltapsi.voila.tsv>;
my @tsv_majiq_hgg  = </Users/naqvia/Desktop/HGG/psi_majiq/*tsv>;
my @tsv_majiq_lgg  = </Users/naqvia/Desktop/LGG/voila_tsv/*non_denovo.lgg_wg_anno.voila_tsv.tsv>;

open(FIL, $samples_file) || die("Cannot Open File $samples_file");
while(<FIL>)
{
  chomp;
  my ($sample_id, $type) = split;
  if($type=~/DIPG/)
  {
    push @dipg_samples, $sample_id;
  }
  elsif($type=~/HGG/)
  {
    push @hgg_samples, $sample_id;
  }
  else{
    push @lgg_samples,$sample_id;
  }
}
close(FIL);

foreach my $file(@tsv_majiq)
{
  my $sample = "";
  if($file=~/(BS\_[\w+\d+]+)\.delta/)
  {
    #print "sample: ",$1,"\t";
    $sample = $1;
  }
  push @samples, $sample;
  open(FIL,$file) || die("Cannot Open File $file");
  while(<FIL>)
  {
      next if ($_=~/^#/);
      my (@content)  = split "\t";
      my $gene_name  = $content[0];
      my $ens        = $content[1];
      my $lsv_coord  = $content[2];
      my $psi_ctrl   = $content[6];
      my $psi_case   = $content[7];
      my $num_exons  = $content[13];
      my $chr        = $content[15];
      my $strand     = $content[16];
      my $exon_coords= $content[18];

      next if($lsv_coord=~/\:t\:/);
      next if($num_exons > 3);
      my ($ctrl_psi, $case_psi);
      my @psi_ctrl = split/;/,$psi_ctrl;
      my @psi_case = split/;/,$psi_case;
      my @exon_coord = split/;/,$exon_coords;

      ## remove complex LSVs
      next unless ( scalar(@psi_ctrl) == (scalar(@exon_coord)-1));

      #print $gene_name,"\t",$lsv_coord,"\t",$psi_ctrl,"\t",$psi_case,"\t";
      #print $exon_coords,"\t";
      #die;

      print $sample,"\t";
      print $ens,"\t",$gene_name,"\t";
      my ($exon_start, $exon_end)       =  split/\-/,$exon_coord[1];
      my ($upstreamES, $upstreamEE)     = split/\-/,$exon_coord[0];
      my ($downstreamES, $downstreamEE) = split/\-/,$exon_coord[2];

      print $lsv_coord,"\t",$gene,"\t",$chr,"\t",$strand,"\t",$exon_start,"\t",$exon_end,"\t";
      print $upstreamEE,"\t",$downstreamES,"\t";
      print $psi_case[0];
      print "\n";

      my $i = 0;
      foreach my $coord(@exon_coord)
      {
        # next unless ( ($exon_coord[$i-1]=~/\d+/) && ($exon_coord[$i+1]=~/\d+/) );
        # print $sample,"\t";
        # print $ens,"\t",$gene_name,"\t";
        # print $lsv_coord,"\t",$gene,"\t",$chr,"\t",$strand,"\t",$coord,"\t";
        # print $exon_coord[$i-1],"\t",$exon_coord[$i-1],"\t";
        # print $psi_case[$i-1];
        # print "\n";
      }
    }
  }
__DATA__
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
