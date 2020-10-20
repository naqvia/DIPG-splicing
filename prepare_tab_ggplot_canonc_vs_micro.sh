#!/bin/sh
## script to prepare tables for ggplot for canonical vs microexon dPSI distributions


grep -f isoform_switch_lsv_ids_list.txt re_analysis/dominant_events_lsvs_simple_total.txt | awk '{if($5>1){ print $0}}' > re_analysis/dominant_events_lsvs_simple_recurrent2.txt
grep -f isoform_switch_lsv_ids_list.txt re_analysis/dominant_events_lsvs_simple_total.txt | awk '{if($5>1){ print $0}}' | awk -F "_" '{print $3}' | awk '{print $1}' | awk -F "-" '{print ($2-$1)+1}' > len.tmp
paste re_analysis/dominant_events_lsvs_simple_recurrent2.txt len.tmp > re_analysis/dominant_events_lsvs_simple_recurrent2.len.txt
cat re_analysis/dominant_events_lsvs_simple_recurrent2.len.txt | awk '{if($8<=30){ print $0 }}' > re_analysis/dominant_events_lsvs_simple_recurrent2.len.micro.txt
cat re_analysis/dominant_events_lsvs_simple_recurrent2.len.txt | awk '{if($8>30){ print $0 }}' > re_analysis/dominant_events_lsvs_simple_recurrent2.len.non_micro.txt
cat re_analysis/dominant_events_lsvs_simple_recurrent2.len.non_micro.txt | awk '{print $1"\t"$2"\tCanonical"}' > re_analysis/dominant_events_lsvs_simple_recurrent2.len.micro_vs_cano.ggplot.txt
cat re_analysis/dominant_events_lsvs_simple_recurrent2.len.micro.txt | awk '{print $1"\t"$2"\tMicroexon"}'>> re_analysis/dominant_events_lsvs_simple_recurrent2.len.micro_vs_cano.ggplot.txt
