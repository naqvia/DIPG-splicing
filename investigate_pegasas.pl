#!/usr/bin/perl

#
# investigate_pegasas.pl
#

my ($signature_gmt, $dir_for_results) = ($ARGV[0],$ARGV[1]);
my @sign_list       = <$dir_for_results/HALLMARK*/*sig_list.txt>;
my @sig_score_files = <$dir_for_results/HALLMARK*/*scores.txt>;

my @groups;
my %gene_count_per_pathway;
my %sig_gene_list_counts;
my %scores_per_pathway;

#print @sign_list;

open(FIL,$signature_gmt) || die("Cannot Open File");
while(<FIL>)
{
  chomp;
  my @cols    = split;
  my $pathway = shift @cols;
  my ($pathway, $remove) = split/\>/,$pathway;
  foreach my $gene(@cols)
  {
    $gene_count_per_pathway{$pathway}++;
  }
}
close(FIL);

foreach my $file(@sign_list)
{
  my $pathway_name = "";
  if($file =~/(\w+)\_sig\_list\.txt/)
  {
    $pathway_name = $1;
  }

  #print "pathway:".$file,"\t*",$pathway_name,"\n";
  push @pathways, $pathway_name;

  open(FIL,$file) || die("Canot Open File $file");
  while(<FIL>)
  {
    chomp;
    $sig_gene_list_counts{$pathway_name}++;
  }
  close(FIL);
}

foreach my $file(@sig_score_files)
{
  my $pathway_name = "";
  if($file =~/(\w+)\.scores.txt/)
  {
    $pathway_name = $1;
  }

  open(FIL,$file) || die("Canot Open File $file");
  while(<FIL>)
  {
    chomp;
    my @cols = split;
    my $sample = $cols[0];
    my $group  = $cols[1];
    my $score  = $cols[2];
    #print "group:".$group," ";

    push @groups,$group;

    push @{$scores_per_pathway{$pathway_name}{$group}}, $score;
    $total_group_count{$group}++;

  }
}

my @groups_uniq = do { my %seen; grep { !$seen{$_}++ } @groups };


foreach my $pathway(@pathways)
{
  print $pathway,"\t";
  my $perc_total = ($sig_gene_list_counts{$pathway}/$gene_count_per_pathway{$pathway}) *100;
  my $rounded =  sprintf("%.0f", $perc_total);
  print $rounded,"\t";
  print $sig_gene_list_counts{$pathway},"\t",$gene_count_per_pathway{$pathway},"\t";

  foreach my $group(@groups_uniq)
  {
    #print "group:".$group,":";

    my $avg_score = &average(\@{$scores_per_pathway{$pathway}{$group}});
  	my $std_score = &stdev(\@{$scores_per_pathway{$pathway}{$group}});

    print $avg_score,"\t",$std_score,"\t";
  }
  print "\n";

}

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
HALLMARK_APOPTOSIS>HALLMARK_APOPTOSIS	CASP3	CASP9	DFFA	CASP7	CFLAR	BIRC3	PMAIP1	CASP8	JUN	BCL2L11	MCL1	IL1B	SPTAN1	DIABLO	BAX	BIK	IL1A	BID	CDKN1A	GADD45A	DDIT3	CDKN1B	TNF	GSN	TNFSF10	CASP6	SQSTM1	FASLG	EGR3	CD44	FAS	IL18	IGFBP6	PRF1	DAP	CCND1	BTG3	F2R	SATB1	BNIP3L	CASP4	TNFRSF12A	CREBBP	RHOB	GPX3	PDGFRB	TSPO	CCND2	XIAP	TIMP1CTNNB1	IRF1	HSPB1	ADD1	TIMP2	BTG2	TIMP3	LEF1	CASP1	GPX1	BCL10	IGF2R	CDC25B	AIFM3	CD38	PPP3R1	HGF	CLU	ATF3	LGALS3	LUM	LMNA	GADD45B	CDK2	IFNB1	RETSAT	SMAD7	SOD1PTK2	ENO2	HMOX1	IER3	BCL2L10	CD2	GCH1	MMP2	VDAC2	TAP1	PLAT	IFNGR1	APP	BRCA1	ROCK1	PSEN1	DCN	PSEN2	SOD2	BMF	EREG	KRT18	TGFB2	RELA	WEE1	RARA	CD14	CD69PEA15	DNAJC3	CASP2	CTH	PLCB2	BMP2	HMGB2	PLPPR4	H1-0	TGFBR3	EBP	TXNIP	ANKH	RHOT2	CYLD	GSTM1	GSR	BGN	BCL2L1	GNA15	MGMT	PPT1	F2	IL6	SC5D	IFITM3	RNASEL	EMP1CAV1	DNM1L	ANXA1	TOP2A	ISG20	SLC20A1	MADD	PPP2R5B	BCAP31	ERBB3	NEDD9	SAT1	PDCD4	BCL2L2	FEZ1	ERBB2	DNAJA1	DAP3	DPYD	NEFH	PAK1	FDXR	GPX4	ETF1	CCNA1	GUCY2D	AVPR1A
