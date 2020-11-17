#!/bin/sh

## generate for ggplot
cat re_analysis/dominant_events_lsvs.total.intersectUnipMod.wo.txt | awk '{print $4"\t"$5"\tModifications"}' > re_analysis/dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat re_analysis/dominant_events_lsvs.total.intersectUnipOther.wo.txt | awk '{print $4"\t"$5"\tOther"}' >> re_analysis/dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat re_analysis/dominant_events_lsvs.total.intersectUnipDisulfBond.wo.txt | awk '{print $4"\t"$5"\tDisulfBond"}' >> re_analysis/dominant_events_lsvs.total.intersectUnip.ggplot.txt
cat re_analysis/dominant_events_lsvs.total.intersectUnipLocSignal.wo.txt | awk '{print $4"\t"$5"\tLocSignal"}' >> re_analysis/dominant_events_lsvs.total.intersectUnip.ggplot.txt
