#!/usr/bin/perl

## make_tab_for_pegasus.pl

my ($samples_file, $expr_file,$pathway_file) = ($ARGV[0], $ARGV[1], $ARGV[2]);
my @pre_dipg_samples;

my (@dipg_nonH3K28_samples,@dipg_H3K28_samples);
my (@lsv);
my (%tab_cols, %psi);
my %already_seen;

## files from MAJIQ output for each cancer Type

my @tsv_majiq_dipg = </Users/naqvia/Desktop/DIPG/psi_analysis/*tsv>;
#my @tsv_majiq_hgg  = </Users/naqvia/Desktop/HGG/psi_majiq/*tsv>;
#my @tsv_majiq_lgg  = </Users/naqvia/Desktop/LGG/voila_tsv/*non_denovo.lgg_wg_anno.voila_tsv.tsv>;

open(FIL, $samples_file) || die("Cannot Open File $samples_file");
while(<FIL>)
{
  chomp;
  my ($sample_id, $type) = split;
  if($type=~/Non/)
  {
    push @dipg_nonH3K28_samples, $sample_id;
  }
  else{
    push @dipg_H3K28_samples,$sample_id;
  }
}
close(FIL);

foreach my $file(@tsv_majiq_dipg)
{
  my $sample = "";
  if($file=~/(BS\_[\w+\d+]+)\.psi/)
  {
    #print "sample: ",$1,"\t";
    $sample = $1;
    $sample=~s/BS_tumor_//;
  }
  push @samples, $sample;
  open(FIL,$file) || die("Cannot Open File $file");
  while(<FIL>)
  {
      #print $_,"\n";
      next if ($_=~/^#/);
      my (@content)  = split "\t";
      my $gene_name  = $content[0];
      my $ens        = $content[1];
      my $lsv_coord  = $content[2];
      #my $psi_ctrl   = $content[6];
      my $psi_case   = $content[3];

      my $num_exons  = $content[10];
      my $chr        = $content[12];
      my $strand     = $content[13];
      my $exon_coords= $content[15];

      next if($lsv_coord=~/\:t\:/);
      next if($num_exons > 3);
      my ($ctrl_psi, $case_psi);
      #my @psi_ctrl = split/;/,$psi_ctrl;
      my @psi_case = split/;/,$psi_case;
      my @exon_coord = split/;/,$exon_coords;

      ## remove complex LSV
      next unless ( scalar(@psi_case) == (scalar(@exon_coord)-1));

      #print $gene_name,"\t",$lsv_coord,"\t",$psi_ctrl,"\t",$psi_case,"\t";
      #print $exon_coords,"\t";
      #die;

      #print "sample:",$sample,"\t";
      #print $ens,"\t",$gene_name,"\t";
      my ($exon_start, $exon_end)       =  split/\-/,$exon_coord[1];
      my ($upstreamES, $upstreamEE)     = split/\-/,$exon_coord[0];
      my ($downstreamES, $downstreamEE) = split/\-/,$exon_coord[2];

      #print $lsv_coord,"\t",$gene,"\t",$chr,"\t",$strand,"\t",$exon_start,"\t",$exon_end,"\t";
      #print $upstreamEE,"\t",$downstreamES,"\t";
      #print $psi_case[0];
      #print "\n";

        if($already_seen{$lsv_coord})
        {
          $psi{$lsv_coord}{$sample} = $psi_case[0];
        }
        else
        {
          #push @{$tab_cols{$lsv_coord}},$ens, $gene_name, $lsv_coord,$chr,$strand,$exon_start,$exon_end,$upstreamEE,$downstreamES;
          push @{$tab_cols{$lsv_coord}},$ens, $gene_name,$chr,$strand,$exon_start,$exon_end,$upstreamEE,$downstreamES;

          $psi{$lsv_coord}{$sample} = $psi_case[0];
          $already_seen{$lsv_coord} = 1;

        }
          push @lsv,$lsv_coord;

          ## AC,GeneName,chr,strand,exonStart,exonEnd,upstreamES,downstreamES,BS_0ZA67BBC,BS_1D6PZNKN BS_1N7MQZGR BS_21ET39G7 BS_2JP7RBMB BS_30VC3R2Q BS_49CJNZ06 BS_5965SFPZ BS_5VPM0F36 BS_68KX6A42 BS_6DCSD5Y6 BS_6M2053M0 BS_7GJT4A6A BS_7WM3MNZ0 BS_8ZY4GST0 BS_97M1E2DW BS_9CA93S6D BS_A7Q8G0Y1 BS_B1C6GZ84 BS_C83TK159 BS_EZ3147MX BS_F0JB4EAK BS_FCDAH728 BS_G3NN392N BS_GWSJ4Z9H BS_H97S5SQN BS_JB43XBCQ BS_JQVAWTTM BS_MKM0EEN1 BS_NB9XXBW6 BS_NEVYM2FP BS_NGSG2KB6 BS_PZVHMSYN BS_Q13FQ8FV BS_QNNX91SM BS_R7NTZR4C BS_R9B92M75 BS_RXN5T5YT BS_W4DB5RP1 BS_W7MFJZ5A BS_WH8G4VFB BS_X4DD4KSZ BS_XGDPK33A BS_XM1AHBDJ BS_YDEVMD24 BS_Z3RCA1T9 BS_ZF6BSFNF
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

  my @lsvs_uniq = do { my %seen; grep { !$seen{$_}++ } @lsv };

print "AC,GeneName,chr,strand,exonStart,exonEnd,upstreamEE,downstreamES";
foreach my $sample(@samples){
  print ",";
  print $sample
}
print "\n";

foreach my $lsv (@lsvs_uniq)
{
  #print "Enter info for LSV:";
  foreach my $item(@{$tab_cols{$lsv}})
  {
    print $item;
    print ",";

  }
  foreach my $sample(@samples){
    if($psi{$lsv}{$sample}){
      my $psi_dec = sprintf("%.2f", $psi{$lsv}{$sample});
      print $psi_dec;
      print ",";

    }
    else{
      print "0";
      print ",";

    }

  }
  print "\n";
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
