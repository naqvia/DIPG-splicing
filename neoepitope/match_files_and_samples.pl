#!usr/bin/perl

my @manifest_files = <*-manifest.csv>;
my @sj_files = <*.SJ.out.tab.gz>;

my (%sample_to_id, %id_to_sample);
my @sj_ids;

foreach my $file($sj_ids)
{
  if($file=~/([\w\d]+\-[\w\d]\-[\w\d]\-[\w\d]\-[\w\d])\.)
  {
    print $file,":",$1;
    push @sj_ids,$1;
  }
}

foreach my $file($manifest_files)
{
  chomp;
  my @cols  = split ",";
  my $name  = $cols[1];
  my $sample_id = "";
  if($_=~/(BS_[\w\d]+)/)
	{
		$sample_id = $1;
	}
  $sample_to_id{$sample_id} = $name;
  $id_to_sample{$name}      = $sample_id;
}

__DATA__
open(FIL,)




__DATA__
id,name,project,Composition,gender,race,ethnicity,Kids First Participant ID,disease_type,sample_id,Tumor Descriptor,sample_type,platform,Kids First Biospecimen ID,primary_site,age_at_diagnosis,aliquot_id,reference_genome,case_id,experimental_strategy
5c7eb4b9e4b0c5cd2e23f935,701f0e2f-49de-4429-81c8-d5162cdc0990.SJ.out.tab.gz,kfdrc-harmonization/sd-bhjxbdqk-11,NULL,Male,Not Reported,Not Reported,PT_RE6AXQM1,Diffuse Astrocytoma,C021_0034,,Tumor,Illumina,BS_6M2053M0,Brain Stem,1460,A19020,GRCh38,P-36,RNAseq
5c7ebb66e4b0c5cd2e240c14,b89862e1-fc40-4779-aaf9-aee0fbc15833.SJ.out.tab.gz,kfdrc-harmonization/sd-bhjxbdqk-11,NULL,Male,White,Hispanic or Latino,PT_NK8A49X5,Diffuse Astrocytoma (WHO gr. 2),C021_0005,,Tumor,Illumina,BS_YDEVMD24,Brain Stem,4745,A08715,GRCh38,P-06,RNAseq
