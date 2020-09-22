#!/usr/bin/perl

#
# clk1_samples.pl
#

my $dir_file = $ARGV[0];
my $clk1_file = $ARGV[1];
#test
my (@non_clk1_samples,@clk1_samples,@splice_id, @nonDMG_samples);
my (%tumor_IJCs, %splice_id_count);

my %clk1_samples_filter;
my %non_clk1_samples_filter;

open(INP,$clk1_file) || die("Cannot Open File");
while(<INP>)
{
  chomp;
  my $sample = $_;
  #print "SAMPLE".$sample."SAMPL\n";
  $clk1_samples_filter{$sample} = $sample;
}
close(INP);

open(FIL,$dir_file) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my $file = $_;
  my $sample = "";
  if($file=~/vs\_(BS\_\w+)\./)
  {
    $sample = $1;
  }
  if($_=~/control\-BR/)
  {
    push @nonDMG_samples,$sample;

  }
  elsif($clk1_samples_filter{$sample})
  {
    push @clk1_samples,$sample;
  }
  else
  {
    push @non_clk1_samples,$sample;
    $non_clk1_samples_filter{$sample} = $sample;
  }

  open(SE_FIL,$file);
  while(<SE_FIL>)
  {
    chomp;
    if($_=~/^\d+/)
    {
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
      #push @samples, $sample;


      push @splice_id, $splice_id;

     # print "*".$sample,"\t",$splice_id,"\t",$tumor_IJC, "\t",$tumor_SJC,"\t",$inc_from_len,"\t",$skip_from_len,"\n";
      $splice_id_count{$splice_id}++;

      ##store values from mis-splicing event
      push @{$spliced_event{$splice_id}{$sample}}, $tumor_IJC, $tumor_SJC,$inc_from_len,$skip_from_len;
      $tumor_IJCs{$splice_id}{$sample} = $tumor_IJC;
      $tumor_SJCs{$splice_id}{$sample} = $tumor_SJC;
      $inc_from_lens{$splice_id} = $inc_from_len;
      $skip_from_lens{$splice_id} = $skip_from_len;
    }
  }
}
  close(SE_FIL);

  my @splice_id_uniq = do { my %seen; grep { !$seen{$_}++ } @splice_id };

  print "SpliceID";
  #print join "\t", @clk1_samples,"***\n";

   foreach my $sample(@non_clk1_samples)
   {
     print ",";
     print $sample;
   }
   foreach my $sample(@clk1_samples)
    {
      print ",";
      print $sample
    }
   print "\n";

  foreach my $splice_id(@splice_id_uniq)
  {
      #next unless $splice_id_count{$splice_id} >=10;
      print $splice_id;

      foreach my $sample(@non_clk1_samples)
      {
        print ",";
        #print $tpm{$sample}{$gene}
        if($tumor_IJCs{$splice_id}{$sample})
        {
            print $tumor_IJCs{$splice_id}{$sample};
        }

        else{
          print "0";
        }
        #my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        #print $rounded;
      }

      foreach my $sample(@clk1_samples)
      {
        print ",";
        if($tumor_IJCs{$splice_id}{$sample})
        {
            print $tumor_IJCs{$splice_id}{$sample};
        }

        else
        {
          print "0";
        }
        #my $rounded = sprintf "%.0f", $tpm{$sample}{$gene};
        #print $rounded;
      }
      print "\n";
    }
