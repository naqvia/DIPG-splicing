#!/usr/bin/perl

##
## make_psi_from_molec_subtype.pl
##

my $dir_file = $ARGV[0];

my (@MT_samples, @WT_samples);
my %tumor_IJCs;
my %tumor_SJCs;
my %splice_id_count;


open(FILM,"dipg_samples_w_molecsubtype.txt") || die("Cannot Open File");
while (<FILM>)
{
  chomp;
  #print $_,"*\n";
  my ($bs_id,$pt_id,$status) = split "\t";
  if($status=~/K28/)
  {
    $status{$bs_id} = "H3K28M";
    #print $_,"\n";
    push @MT_samples, $bs_id;
  }
  else{
    $status{$bs_id} = "WT";
    push @WT_samples, $bs_id;
  }
}
close(FILM);

open(FIL,$dir_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my $file = $_;
  my $sample = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    #print "samp*".$1,"*\n";
    $sample = $1;
  }
  open(SE_FIL,$file);
  while(<SE_FIL>)
  {
    chomp;
    if($_=~/^\d+/)
    {
      #print "*".$_,"\n";

      my @cols = split "\t";
      my $gene = $cols[2];
      my $chr  = $cols[3];
      my $str  =   $cols[4];

      my $exon_start    = $cols[5];
      my $exon_end      = $cols[6];
      my $tumor_IJC     = $cols[14];
      my $tumor_SJC     = $cols[15];
      my $inc_from_len  = $cols[16];
      my $skip_from_len = $cols[17];
      my $pval = $cols[18];
      my $thr_diff = $cols[-1];

      ## only look at strong changes
      next unless (abs($thr_diff) >= .10);
      next unless ($pval<=0.05);

      next unless ($tumor_IJC >=10);
      next unless ($tumor_SJC >=10);


      #print "*",$_,"\n";
      my $splice_id = $gene."_".$chr.":".$exon_start."-".$exon_end;
      $splice_id=~s/\"//g;
      next unless ($pval<=0.05);
      push @samples, $sample;
      push @splice_id, $splice_id;
      $samples_filter{$sample} = $sample;

    #  print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";

      push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      $splice_id_count{$splice_id}++;
      #$spliced_event{$splice_id}{$sample}
      $tumor_IJCs{$splice_id}{$sample} = $tumor_IJC;
      $tumor_SJCs{$splice_id}{$sample} = $tumor_SJC;
      $inc_from_lens{$splice_id} = $inc_from_len;
      $skip_from_lens{$splice_id} = $skip_from_len;
    }
  }
}
  close(SE_FIL);

  my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_id };
  # ## $tumor_IJCs{$splice_id}{$sample}
  # print "SpliceID";
  # foreach my $sample(@WT_samples)
  # {
  #   print ",";
  #   print $sample
  # }
  # foreach my $sample(@MT_samples)
  # {
  #   print ",";
  #   print $sample;
  # }
  # print "\n";

  foreach my $splice_id(@splice_id_uniq)
  {
      next unless $splice_id_count{$splice_id} >=10;
      #print $splice_id;

      # foreach my $sample(@WT_samples)
      # {
      #   print ",";
      #   #print $tpm{$sample}{$gene}
      #   if($tumor_IJCs{$splice_id}{$sample})
      #   {
      #       print $tumor_IJCs{$splice_id}{$sample};
      #   }
      #
      #   else{
      #     #print "0";
      #   }
      #   #my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      #   #print $rounded;
      # }
      #
      # foreach my $sample(@MT_samples)
      # {
      #   print ",";
      #   if($tumor_IJCs{$splice_id}{$sample})
      #   {
      #       print $tumor_IJCs{$splice_id}{$sample};
      #   }
      #
      #   else
      #   {
      #     #print "0";
      #   }
      #   #my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
      #   #print $rounded;
      # }

      foreach my $sample(@WT_samples)
      {
        #print $tpm{$sample}{$gene}
        if($tumor_SJCs{$splice_id}{$sample})
        {
            #print $tumor_IJCs{$splice_id}{$sample};
            print $splice_id;
            print "\tWT\n";
        }

      }

      foreach my $sample(@MT_samples)
      {
        if($tumor_SJCs{$splice_id}{$sample})
        {
          #  print $tumor_IJCs{$splice_id}{$sample};
          print $splice_id;
          print "\tMT\n";
        }

      }

    }
