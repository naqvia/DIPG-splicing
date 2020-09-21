#!/bin/sh

#perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | awk '{print $6"\t"$1"\t"$2"\t"$7}'| perl -pe 's/(chr[\d+XY]+)\:/$1\t/g' | awk -F "-" '{print $1"\t"$2"-"$3"-"$4}' | awk '{if($6 == "+"){ print $0 } else {print $0"-" }}' | perl -pe 's/\*//' > output/dominant_changes_lsv.bed
#perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | perl -pe 's/\*//'  > output/dominant_changes_lsv.txt
#perl scripts/extract_recurrent_dominant_lsvs.pl | grep "*" | perl scripts/reformat_for_bed.pl >  output/dominant_changes_lsv.bed

## Unipro cross reference
#bedtools intersect -wo -a output/dominant_changes_lsv.bed -b unipMod.hg38.col.bed > output/dominant_events_lsvs.total.intersectUnipMod.wo.txt
#bedtools intersect -wo -a output/dominant_changes_lsv.bed -b unipOther.hg38.col.bed > output/dominant_events_lsvs.total.intersectUnipOther.wo.txt
#bedtools intersect -wo -a output/dominant_changes_lsv.bed -b unipDisulfBond.sorted.col.bed > output/dominant_events_lsvs.total.intersectUnipDisulfBond.wo.txt
#bedtools intersect -wo -a output/dominant_changes_lsv.bed -b unipLocSignal.sorted.col.bed > output/dominant_events_lsvs.total.intersectUnipLocSignal.wo.txt
#bedtools intersect -wo -a output/dominant_changes_lsv.bed -b unipDomain.hg38.col.bed > output/dominant_events_lsvs.total.intersectUnipDomain.wo.txt

## generate for ggplot
cat output/dominant_events_lsvs.total.intersectUnipMod.wo.txt | awk '{print $4"\t"$5"\tModifications"}' > output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat output/dominant_events_lsvs.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat output/dominant_events_lsvs.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat output/dominant_events_lsvs.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt

# cat output/dominant_events_lsvs.total.intersectUnipDomain.wo.txt | awk '{print $4"\t"$5"\tDomain"}' >> output/dominant_events_lsvs.total.intersectUnip.ggplot.txt
