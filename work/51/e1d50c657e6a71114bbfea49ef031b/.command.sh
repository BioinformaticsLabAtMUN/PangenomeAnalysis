#!/bin/bash -ue
python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/heaps_analysis.py \
    --matrix Strep_cdhit_strain_by_gene.npz \
    --labels Strep_cdhit_strain_by_gene.npz.labels.txt \
    --output . \
    --iterations 100

mv heaps_law_results.csv   Strep_heaps_law_results.csv
mv heaps_law_iterations.csv Strep_heaps_law_iterations.csv
mv heaps_law_data.npz      Strep_heaps_law_data.npz
mv heaps_law_plot.png      Strep_heaps_law_plot.png
