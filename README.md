## DIPG_aberrant-splicing  

1. General overview
2. Splicing factors dysregulation
 -SRSF family and poison exon mis-splicing
 -Expression analysis and methylation results for splicing factors
3. Cross reference with Uniprot domain information and functional sites
 -identify loss/gain of domains
 -loss/gain of functional sites
4. "isoform switching" events
5. Histone and chromatin regulation
 -General overview of their aberrant splicing (summary type plots)
 -Focusing on tumor suppressor SMARCA4 and DPF2 (since these are prevalent, recurrent and a part of the BAF complex)
 -Establish that their targets are also compromised
 -Strong isoform-switching splicing events on tsg and oncogenes
  -general plot with all, and then focusing on recurrent ones
  -select few tumor suppressors and oncogenes (those that could have functional consequences)
6. Micro-exons
 -summary plot
 -recurrent ones, that are isoform switching, and that have associated functional sites affected
 -characterize out-of-frame micro-exons and focus on oncogenic BAK1
 -expand to de-novo micro-exons

## pipeline flow for protein domain and functional sites
##run script to find dominant isoform and generate bed files for exons
perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | awk '{print $6"\t"$1"\t"$2"\t"$7}'| perl -pe 's/(chr[\d+XY]+)\:/$1\t/g' | awk -F "-" '{print $1"\t"$2"-"$3"-"$4}' | awk '{if($6 == "+"){ print $0 } else {print $0"-" }}' | perl -pe 's/\*//' > dominant_events_lsvs.total.bed


## pre-process unipro bed files (downloaded from UCSC):
sort -k1,1 -k2,2n ~/Desktop/reference_data/pfam.hg38.bed > ~/Desktop/reference_data/pfam.hg38.sorted.bed
bedtools merge -s -c 4,4,6 -header -o distinct,count_distinct,distinct -i ~/Desktop/reference_data/pfam.hg38.sorted.bed > pfam.hg38.col.bed

## Unipro cross reference
bedtools intersect -wo -a dominant_events_lsvs.total.bed -b unipMod.hg38.col.bed > dominant_events_lsvs.total.intersectUnipMod.wo.txt
bedtools intersect -wo -a dominant_events_lsvs.total.bed -b unipOther.hg38.col.bed > dominant_events_lsvs.total.intersectUnipOther.wo.txt
bedtools intersect -wo -a dominant_events_lsvs.total.bed -b unipDisulfBond.sorted.col.bed > dominant_events_lsvs.total.intersectUnipDisulfBond.wo.txt
bedtools intersect -wo -a dominant_events_lsvs.total.bed -b unipLocSignal.sorted.col.bed > dominant_events_lsvs.total.intersectUnipLocSignal.wo.txt
bedtools intersect -wo -a dominant_events_lsvs.total.bed -b unipDomain.hg38.col.bed > dominant_events_lsvs.total.intersectUnipDomain.wo.txt

## generate for ggplot
  cat dominant_events_lsvs.100perc.intersectUnipMod.wo.txt | awk '{print $4"\t"$5"\tModifications"}' > dominant_events_lsvs.100perc.intersectUnip.ggplot.txt
  545  cat dominant_events_lsvs.100perc.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' >> dominant_events_lsvs.100perc.intersectUnip.ggplot.txt
  546  cat dominant_events_lsvs.100perc.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' >> dominant_events_lsvs.100perc.intersectUnip.ggplot.txt
  547  cat dominant_events_lsvs.100perc.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' >> dominant_events_lsvs.100perc.intersectUnip.ggplot.txt

## neoepitope project
1. perl ../../../scripts/extract_juncs_coords_bed.pl
2. perl ../scripts/neoepitope/reformat_bed_for_junc.perl
3. perl ../scripts/neoepitope/make_combined_tab.pl > DIPG_psi_neo_list_021720.txt
4. perl ../scripts/neoepitope/make_final_tab_for_neojunc.pl
