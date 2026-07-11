#!/bin/bash -ue
echo "Starting core genome analysis..."
echo "Input matrix: Strep_cdhit_strain_by_gene.npz"
echo "Input labels: Strep_cdhit_strain_by_gene.npz.labels.txt"

python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/core_genome_analysis.py \
    --matrix Strep_cdhit_strain_by_gene.npz \
    --labels Strep_cdhit_strain_by_gene.npz.labels.txt \
    --output-prefix Strep

echo "Core genome analysis completed"
