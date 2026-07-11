#!/bin/bash -ue
echo "=== Gene-Level Structural Analysis for cdhit ==="

export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen

python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/analyze_gene_structure.py \
    --gene-matrix Strep_cdhit_strain_by_gene.npz \
    --gene-labels Strep_cdhit_strain_by_gene.npz.labels.txt \
    --core-genes Strep_core_genes.txt \
    --metadata my_proteome_metadata.tsv \
    --annotations Strep_cdhit_core_genes_annotated.tsv \
    --output-dir . \
    --max-strains-plot 50

# Rename all output files to include the method prefix
for file in *.tsv *.txt; do
    if [ -f "$file" ] && [[ "$file" != "Strep_cdhit_"* ]]; then
        base_name=$(basename "$file" | sed 's/\.[^.]*$//')
        extension=$(basename "$file" | sed 's/.*\.//')
        mv "$file" "Strep_cdhit_${base_name}.${extension}"
    fi
done

for png_file in *.png; do
    if [ -f "$png_file" ] && [[ "$png_file" != "Strep_cdhit_"* ]]; then
        base_name=$(basename "$png_file" .png)
        mv "$png_file" "Strep_cdhit_${base_name}.png"
    fi
done

for html_file in *.html; do
    if [ -f "$html_file" ] && [[ "$html_file" != "Strep_cdhit_"* ]]; then
        base_name=$(basename "$html_file" .html)
        mv "$html_file" "Strep_cdhit_${base_name}.html"
    fi
done

echo "=== Gene Structural Analysis Complete ==="
