#!/bin/sh

cat re_analysis/dominant_events_lsvs.total.intersectUnipMod.wo.txt   | awk '{print $4"\t"$5"\tModifications"}' > dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat re_analysis/dominant_events_lsvs.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' >> dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat re_analysis/dominant_events_lsvs.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' >> dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat re_analysis/dominant_events_lsvs.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' >> dominant_events_lsvs.total.intersectUnip.ggplot.txt
