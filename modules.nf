// =============================================================================
//  modules.nf  —  all shared processes
//  1. $projectDir is a Nextflow built-in that always equals the directory
//     containing main.nf, regardless of who runs the pipeline or from where.
//     params.baseDir was set to $PWD (wherever the user *ran* the pipeline),
//     which means it would break the moment someone ran it from a different
//     directory, or handed it to another person.
//
//  2. Each process now has a   conda   directive at the top.
//     Nextflow reads that directive, builds/caches the environment once,
//     and activates it automatically before the script runs.
//     This means:
//       - No custom shell script needed
//       - No assumption about where conda is installed
//       - Works on any machine, cluster, or cloud that has conda
//
//  3. Nextflow manages activation and deactivation itself.
//
//  4. Pckages are now declared in envs/pangenome_env.yml.
// =============================================================================


// -----------------------------------------------------------------------------
//  renameSequences
// -----------------------------------------------------------------------------
process renameSequences {
   
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.renamed_directory}/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/renamed/${method}"

    input:
    tuple path(cluster_file), path(fasta_file)
    val method
    path shared_headers_file
    path mapping_file

    output:
    tuple path("${params.name_prefix}_${method}_renamed.fasta"),
          path("${params.name_prefix}_${method}_allele_names.tsv"),
          val(method)

    script:
 
    """
    SCRIPT_PATH="$projectDir/scripts/unified_rename_sequences_simple_representative.py"

    echo -e "UniProt_ID\\tPDB_file\\tStatus" > dummy_mapping.tsv

    if [[ "${method}" != "foldseek" ]]; then
        MAPPING_FILE=dummy_mapping.tsv
    else
        MAPPING_FILE=${mapping_file}
    fi

    python \$SCRIPT_PATH \\
        ${cluster_file} \\
        ${fasta_file} \\
        ${params.name_prefix}_${method}_renamed.fasta \\
        ${params.name_prefix}_${method}_allele_names.tsv \\
        ${params.name_prefix} \\
        ${method} \\
        ${shared_headers_file} \\
        \$MAPPING_FILE
    """
}


// -----------------------------------------------------------------------------
//  generatePangenomeTables
// -----------------------------------------------------------------------------
process generatePangenomeTables {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/pangenome_tables/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/pangenome_tables/${method}"

    input:
    tuple path(renamed_fasta), path(allele_names), val(method)
    path input_directory

    output:
    tuple path("${params.name_prefix}_${method}_strain_by_allele.npz"),
          path("${params.name_prefix}_${method}_strain_by_gene.npz"),
          path("${params.name_prefix}_${method}_strain_by_allele.npz.labels.txt"),
          path("${params.name_prefix}_${method}_strain_by_gene.npz.labels.txt"),
          val(method)

    script:
    """
    mkdir -p input_fasta

    find ${input_directory} -follow -type f \\( -iname "*.fa" -o -iname "*.faa" -o -iname "*.fasta" -o -iname "*.fas" \\) 2>/dev/null | \\
      while read file; do
        base=\$(basename "\$file")
        if [[ "\$base" != *"renamed"* && "\$base" != *"${method}"* ]]; then
            ln -sf "\$(readlink -f "\$file")" input_fasta/
        fi
      done

    python $projectDir/scripts/build_pangenome_tables.py \\
        --input-dir input_fasta \\
        --output-dir . \\
        --name ${params.name_prefix}_${method} \\
        --allele-names ${allele_names} \\
        --shared-headers ${params.shared_headers_file}
    """
}


// -----------------------------------------------------------------------------
//  validatePangenome
// -----------------------------------------------------------------------------
process validatePangenome {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/validation/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/validate_pangenome/${method}"
    debug true

    input:
    tuple path(allele_matrix), path(gene_matrix),
          path(allele_labels), path(gene_labels), val(method)
    tuple path(renamed_faa), path(allele_names)
    path input_directory

    output:
    path "${params.name_prefix}_${method}_validation_summary.txt"

    script:
    """
    export PYTHONUNBUFFERED=1

    python $projectDir/scripts/validate_pangenome.py \\
        --gene-matrix ${gene_matrix} \\
        --allele-matrix ${allele_matrix} \\
        --input-dir ${input_directory} \\
        --allele-names ${allele_names} \\
        --workers 1 \\
        --batch-size 50 \\
        --output-dir . 2>&1 | tee ${params.name_prefix}_${method}_validation_summary.txt
    """
}


// -----------------------------------------------------------------------------
//  viewPangenome
// -----------------------------------------------------------------------------
process viewPangenome {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/pangenome_visualizations/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/visualize_pangenome/${method}"

    input:
    tuple path(allele_matrix), path(gene_matrix), path(allele_labels), path(gene_labels), val(method)

    output:
    path "pangenome_viz.log"
    path "*_gene_*.png" optional true
    path "*_no_heatmap.txt" optional true

    script:
    """
    DATA_SUMMARY_PARAM=""
    if [ -f "${params.data_summary_file}" ]; then
        echo "Using data_summary.tsv for organism mapping"
        DATA_SUMMARY_PARAM="--data-summary ${params.data_summary_file}"
    fi

    python $projectDir/scripts/pan_viz.py \\
        --input-dir . \\
        --output-dir . \\
        --name ${params.name_prefix}_${method} \\
        --max-features 150 \\
        --core-boost 4 \\
        --figure-width 18 \\
        --figure-height 14 \\
        \$DATA_SUMMARY_PARAM 2>&1 | tee pangenome_viz.log

    if [ \$? -ne 0 ]; then
        echo "Visualization failed, creating fallback file"
        touch "${params.name_prefix}_${method}_no_heatmap.txt"
    fi
    """
}


// -----------------------------------------------------------------------------
//  analyzeHeapsLaw
// -----------------------------------------------------------------------------
process analyzeHeapsLaw {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/heaps_analysis/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/heaps_law/${method}"

    input:
    tuple path(gene_matrix), path(gene_labels)
    val method

    output:
    tuple path("${params.name_prefix}_heaps_law_results.csv"),
          path("${params.name_prefix}_heaps_law_iterations.csv"),
          path("${params.name_prefix}_heaps_law_data.npz"),
          path("${params.name_prefix}_heaps_law_plot.png")

    script:
    """
    python $projectDir/scripts/heaps_analysis.py \\
        --matrix ${gene_matrix} \\
        --labels ${gene_labels} \\
        --output . \\
        --iterations 100

    mv heaps_law_results.csv   ${params.name_prefix}_heaps_law_results.csv
    mv heaps_law_iterations.csv ${params.name_prefix}_heaps_law_iterations.csv
    mv heaps_law_data.npz      ${params.name_prefix}_heaps_law_data.npz
    mv heaps_law_plot.png      ${params.name_prefix}_heaps_law_plot.png
    """
}


// -----------------------------------------------------------------------------
//  analyzeAndValidateCoreGenome
// -----------------------------------------------------------------------------
process analyzeAndValidateCoreGenome {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/core_genome/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/core/${method}"

    input:
    tuple path(gene_matrix), path(gene_labels)
    val method

    output:
    tuple path("${params.name_prefix}_core_genes.txt"),
          path("${params.name_prefix}_core_matrix.npz"),
          path("${params.name_prefix}_frequency_estimates.csv"),
          path("${params.name_prefix}_beta_binomial_results.csv"),
          path("${params.name_prefix}_beta_binomial_fit.png", optional: true)

    script:
    """
    echo "Starting core genome analysis..."
    echo "Input matrix: ${gene_matrix}"
    echo "Input labels: ${gene_labels}"

    python $projectDir/scripts/core_genome_analysis.py \\
        --matrix ${gene_matrix} \\
        --labels ${gene_labels} \\
        --output-prefix ${params.name_prefix}

    echo "Core genome analysis completed"
    """
}


// -----------------------------------------------------------------------------
//  extractDominantAllelesWithCore
// -----------------------------------------------------------------------------
process extractDominantAllelesWithCore {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/dominant_alleles/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/dominant_alleles/${method}"

    input:
    tuple path(allele_matrix), path(gene_matrix),
          path(allele_labels), path(gene_labels), val(method),
          path(renamed_fasta), path(allele_names_tsv), path(core_genes_list)

    output:
    tuple path("${params.name_prefix}_${method}_all_dominant_alleles.faa"),
          path("${params.name_prefix}_${method}_all_dominant_summary.tsv"),
          path("${params.name_prefix}_${method}_core_dominant_alleles.faa"),
          path("${params.name_prefix}_${method}_core_dominant_summary.tsv"),
          val(method)

    script: 
    """
    python $projectDir/scripts/extract_dominant_alleles.py \\
        --allele-matrix ${allele_matrix} \\
        --allele-labels ${allele_labels} \\
        --allele-names ${allele_names_tsv} \\
        --allele-faa ${renamed_fasta} \\
        --core-genes ${core_genes_list} \\
        --all-dominant-faa ${params.name_prefix}_${method}_all_dominant_alleles.faa \\
        --all-dominant-summary ${params.name_prefix}_${method}_all_dominant_summary.tsv \\
        --core-dominant-faa ${params.name_prefix}_${method}_core_dominant_alleles.faa \\
        --core-dominant-summary ${params.name_prefix}_${method}_core_dominant_summary.tsv \\
        --max-alternatives ${params.max_alternatives_per_cluster}
    """
}


// -----------------------------------------------------------------------------
//  annotateSmartDominantAlleles
// -----------------------------------------------------------------------------
process annotateSmartDominantAlleles {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/annotations/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/annotations/${method}"

    time '4h'
    memory '8 GB'

    input:
    tuple path(all_summary), path(core_summary), val(method)

    output:
    path "${params.name_prefix}_${method}_core_genes_annotated.tsv",      emit: core_annotations
    path "${params.name_prefix}_${method}_accessory_genes_annotated.tsv", emit: accessory_annotations, optional: true

    script:
    if (params.database == "ncbi")
    """
    echo "Using NCBI annotation with local GTFs"
    python $projectDir/scripts/annotate_ncbi_with_gtf.py \\
        --input-summary ${all_summary} \\
        --core-summary ${core_summary} \\
        --gtf-directory ${params.gtf_directory} \\
        --annotate-scope ${params.annotate_scope} \\
        --output-core-annotations ${params.name_prefix}_${method}_core_genes_annotated.tsv \\
        --output-core-merged ${params.name_prefix}_${method}_core_genes_with_annotations.tsv \\
        --output-accessory-annotations ${params.name_prefix}_${method}_accessory_genes_annotated.tsv \\
        --output-accessory-merged ${params.name_prefix}_${method}_accessory_genes_with_annotations.tsv
    """

    else if (params.database == "uniprot")
    """
    echo "Using UniProt batch annotation"

    python $projectDir/scripts/annotate_uniprot_batch.py \\
        --input-summary ${all_summary} \\
        --core-summary ${core_summary} \\
        --annotate-scope ${params.annotate_scope} \\
        --output-core-annotations ${params.name_prefix}_${method}_core_genes_annotated.tsv \\
        --output-core-merged ${params.name_prefix}_${method}_core_genes_with_annotations.tsv \\
        --output-accessory-annotations ${params.name_prefix}_${method}_accessory_genes_annotated.tsv \\
        --output-accessory-merged ${params.name_prefix}_${method}_accessory_genes_with_annotations.tsv \\
        --max-alternatives ${params.max_alternatives_per_cluster}
    """

    else
    error "Invalid database parameter: ${params.database}. Must be 'ncbi' or 'uniprot'"
}


// -----------------------------------------------------------------------------
//  analyzeGeneStructure
// -----------------------------------------------------------------------------
process analyzeGeneStructure {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/gene_structure_analysis/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/gene_structure_analysis/${method}"
    debug true

    input:
    tuple path(gene_matrix_npz), path(gene_labels), path(core_genes),
          path(proteome_metadata_file), path(core_annotations), val(method)

    output:
    tuple path("${params.name_prefix}_${method}_core_genome_structural_analysis.tsv"),
          path("${params.name_prefix}_${method}_unique_genes_per_strain_analysis.tsv"),
          path("${params.name_prefix}_${method}_gene_frequency_analysis.tsv"),
          path("${params.name_prefix}_${method}_missing_core_genes_structural.tsv"),
          path("${params.name_prefix}_${method}_frequently_missing_core_genes_structural.tsv"),
          path("${params.name_prefix}_${method}_structural_analysis_summary.txt"),
          path("${params.name_prefix}_${method}_interactive_core_coverage_with_outlier_names.html", optional: true),
          path("${params.name_prefix}_${method}_color_coded_full_species_names.png"),
          path("${params.name_prefix}_${method}_pangenome_composition_simplified.png"),
          path("${params.name_prefix}_${method}_all_strains_missing_core_genes.png"),
          path("${params.name_prefix}_${method}_most_frequently_missing_core_genes_detailed.png"),
          path("${params.name_prefix}_${method}_missing_core_genes_summary_table.tsv"),
          path("${params.name_prefix}_${method}_missing_core_genes_summary_table.html"),
          path("${params.name_prefix}_${method}_unique_genes_distribution.png"),
          val(method)

    script:
    """
    echo "=== Gene-Level Structural Analysis for ${method} ==="

    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen

    python $projectDir/scripts/analyze_gene_structure.py \\
        --gene-matrix ${gene_matrix_npz} \\
        --gene-labels ${gene_labels} \\
        --core-genes ${core_genes} \\
        --metadata ${proteome_metadata_file} \\
        --annotations ${core_annotations} \\
        --output-dir . \\
        --max-strains-plot ${params.max_strains_missing_plot}

    # Rename all output files to include the method prefix
    for file in *.tsv *.txt; do
        if [ -f "\$file" ] && [[ "\$file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$file" | sed 's/\\.[^.]*\$//')
            extension=\$(basename "\$file" | sed 's/.*\\.//')
            mv "\$file" "${params.name_prefix}_${method}_\${base_name}.\${extension}"
        fi
    done

    for png_file in *.png; do
        if [ -f "\$png_file" ] && [[ "\$png_file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$png_file" .png)
            mv "\$png_file" "${params.name_prefix}_${method}_\${base_name}.png"
        fi
    done

    for html_file in *.html; do
        if [ -f "\$html_file" ] && [[ "\$html_file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$html_file" .html)
            mv "\$html_file" "${params.name_prefix}_${method}_\${base_name}.html"
        fi
    done

    echo "=== Gene Structural Analysis Complete ==="
    """
}


// -----------------------------------------------------------------------------
//  analyzeFunctionalCore
// -----------------------------------------------------------------------------
process analyzeFunctionalCore {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/functional_core_analysis/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/functional_core_analysis/${method}"
    debug true

    input:
    tuple path(gene_matrix_npz), path(gene_labels), path(core_genes),
          path(core_annotations), path(proteome_metadata_file), val(method)

    output:
    tuple path("${params.name_prefix}_${method}_core_go_biological_process.tsv"),
          path("${params.name_prefix}_${method}_core_go_molecular_function.tsv"),
          path("${params.name_prefix}_${method}_core_protein_families.tsv"),
          path("${params.name_prefix}_${method}_strain_missing_functional_analysis.tsv"),
          path("${params.name_prefix}_${method}_frequently_missing_go_bp.tsv"),
          path("${params.name_prefix}_${method}_frequently_missing_protein_families.tsv"),
          path("${params.name_prefix}_${method}_*_genes_go_bp.tsv"),
          path("${params.name_prefix}_${method}_*_genes_protein_families.tsv"),
          path("${params.name_prefix}_${method}_functional_similarity_analysis.tsv", optional: true),
          path("${params.name_prefix}_${method}_functional_analysis_summary.txt"),
          path("*.png"),
          val(method)

    script:
    """
    echo "=== Functional Core Genome Analysis for ${method} ==="

    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen

    python $projectDir/scripts/analyze_functional_core.py \\
        --gene-matrix ${gene_matrix_npz} \\
        --gene-labels ${gene_labels} \\
        --core-genes ${core_genes} \\
        --annotations ${core_annotations} \\
        --metadata ${proteome_metadata_file} \\
        --output-dir .

    for file in *.tsv *.txt; do
        if [ -f "\$file" ] && [[ "\$file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$file" | sed 's/\\.[^.]*\$//')
            extension=\$(basename "\$file" | sed 's/.*\\.//')
            mv "\$file" "${params.name_prefix}_${method}_\${base_name}.\${extension}"
        fi
    done

    for png_file in *.png; do
        if [ -f "\$png_file" ] && [[ "\$png_file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$png_file" .png)
            mv "\$png_file" "${params.name_prefix}_${method}_\${base_name}.png"
        fi
    done

    echo "=== Functional Core Analysis Complete ==="
    """
}


// -----------------------------------------------------------------------------
//  clusterGOTermsWithGOATools
// -----------------------------------------------------------------------------
process clusterGOTermsWithGOATools {
    // This process needs GOATools which has different dependencies, so it keeps
    // its own separate environment file.
    conda "$projectDir/envs/goatools_env.yml"

    publishDir "${params.enhanced_output}/goatools_clustering/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/goatools_clustering/${method}"

    input:
    tuple path(core_annotations), path(accessory_annotations), val(method)

    output:
    path "${params.name_prefix}_${method}_goatools_bp_clusters.tsv",          emit: bp_clusters
    path "${params.name_prefix}_${method}_goatools_mf_clusters.tsv",          emit: mf_clusters
    path "${params.name_prefix}_${method}_goatools_summary.txt",              emit: summary
    path "*.png",                                                              emit: plots, optional: true
    path "${params.name_prefix}_${method}_core_goatools_bp_clusters.tsv",     emit: core_bp, optional: true
    path "${params.name_prefix}_${method}_core_goatools_mf_clusters.tsv",     emit: core_mf, optional: true
    path "${params.name_prefix}_${method}_accessory_goatools_bp_clusters.tsv", emit: acc_bp, optional: true
    path "${params.name_prefix}_${method}_accessory_goatools_mf_clusters.tsv", emit: acc_mf, optional: true

    script:
    """
    echo "=== GO Term Clustering with GOATools (Resnik similarity) for ${method} ==="

    python $projectDir/scripts/cluster_go_terms_goatools.py \\
        --core-annotations ${core_annotations} \\
        --accessory-annotations ${accessory_annotations} \\
        --output-prefix ${params.name_prefix}_${method} \\
        --method ${method} \\
        --threshold 0.7

    echo "=== GOATools Clustering Complete ==="
    ls -la *goatools*
    """
}


// -----------------------------------------------------------------------------
//  analyzeRevigoFunctionalCore
// -----------------------------------------------------------------------------
process analyzeRevigoFunctionalCore {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/revigo_functional_core_analysis/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/revigo_functional_core_analysis/${method}"
    debug true

    input:
    tuple path(gene_matrix_npz), path(gene_labels), path(core_genes),
          path(core_annotations), path(proteome_metadata_file),
          path(revigo_bp_file), path(revigo_mf_file), val(method)

    output:
    path "${params.name_prefix}_${method}_revigo_*.tsv"
    path "${params.name_prefix}_${method}_revigo_*.txt"
    path "${params.name_prefix}_${method}_strain_*.html"
    path "${params.name_prefix}_${method}_strain_*.tsv"
    path "*.png"

    script:
    """
    echo "=== Revigo-Based Functional Core Analysis for ${method} ==="

    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen

    python $projectDir/scripts/analyze_revigo_functional_core.py \\
        --gene-matrix ${gene_matrix_npz} \\
        --gene-labels ${gene_labels} \\
        --core-genes ${core_genes} \\
        --annotations ${core_annotations} \\
        --metadata ${proteome_metadata_file} \\
        --revigo-bp ${revigo_bp_file} \\
        --revigo-mf ${revigo_mf_file} \\
        --output-dir . \\
        --debug

    for file in *.tsv *.txt *.html; do
        if [ -f "\$file" ] && [[ "\$file" != "${params.name_prefix}_${method}_"* ]]; then
            if [[ "\$file" == "revigo_debug_"* ]]; then
                base_name=\$(echo "\$file" | sed 's/revigo_debug_//')
                mv "\$file" "${params.name_prefix}_${method}_revigo_debug_\${base_name}"
            elif [[ "\$file" == "strain_"* ]]; then
                mv "\$file" "${params.name_prefix}_${method}_\$file"
            elif [[ "\$file" == "revigo_"* ]]; then
                mv "\$file" "${params.name_prefix}_${method}_\$file"
            else
                base_name=\$(basename "\$file" | sed 's/\\.[^.]*\$//')
                extension=\$(basename "\$file" | sed 's/.*\\.//')
                mv "\$file" "${params.name_prefix}_${method}_\${base_name}.\${extension}"
            fi
        fi
    done

    for png_file in *.png; do
        if [ -f "\$png_file" ] && [[ "\$png_file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$png_file" .png)
            mv "\$png_file" "${params.name_prefix}_${method}_\${base_name}.png"
        fi
    done

    echo "=== Analysis Complete ==="
    ls -la ${params.name_prefix}_${method}_* *.png
    """
}


// -----------------------------------------------------------------------------
//  analyzeClusteredFunctionalCore
// -----------------------------------------------------------------------------
process analyzeClusteredFunctionalCore {
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/clustered_functional_core_analysis/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/clustered_functional_core_analysis/${method}"
    debug true

    input:
    tuple path(gene_matrix_npz), path(gene_labels), path(core_genes),
          path(core_annotations), path(proteome_metadata_file),
          path(goatools_bp_file), path(goatools_mf_file), path(allele_names_file), val(method)

    output:
    path "${params.name_prefix}_${method}_clustered_*.tsv"
    path "${params.name_prefix}_${method}_clustered_*.txt"
    path "${params.name_prefix}_${method}_strain_*.html"
    path "${params.name_prefix}_${method}_strain_*.tsv"
    path "*.png"

    script:
    """
    echo "=== GOATools Clustering-Based Functional Core Analysis for ${method} ==="

    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen

    python $projectDir/scripts/clustered_functional_core.py \\
        --gene-matrix ${gene_matrix_npz} \\
        --gene-labels ${gene_labels} \\
        --core-genes ${core_genes} \\
        --annotations ${core_annotations} \\
        --metadata ${proteome_metadata_file} \\
        --goatools-bp ${goatools_bp_file} \\
        --goatools-mf ${goatools_mf_file} \\
        --allele-names ${allele_names_file} \\
        --output-dir . \\
        --debug

    for file in *.tsv *.txt *.html; do
        if [ -f "\$file" ] && [[ "\$file" != "${params.name_prefix}_${method}_"* ]]; then
            if [[ "\$file" == "revigo_debug_"* ]]; then
                base_name=\$(echo "\$file" | sed 's/revigo_debug_//')
                mv "\$file" "${params.name_prefix}_${method}_revigo_debug_\${base_name}"
            elif [[ "\$file" == "strain_"* ]]; then
                mv "\$file" "${params.name_prefix}_${method}_\$file"
            elif [[ "\$file" == "revigo_"* ]]; then
                mv "\$file" "${params.name_prefix}_${method}_\$file"
            else
                base_name=\$(basename "\$file" | sed 's/\\.[^.]*\$//')
                extension=\$(basename "\$file" | sed 's/.*\\.//')
                mv "\$file" "${params.name_prefix}_${method}_\${base_name}.\${extension}"
            fi
        fi
    done

    for png_file in *.png; do
        if [ -f "\$png_file" ] && [[ "\$png_file" != "${params.name_prefix}_${method}_"* ]]; then
            base_name=\$(basename "\$png_file" .png)
            mv "\$png_file" "${params.name_prefix}_${method}_\${base_name}.png"
        fi
    done

    echo "=== Analysis Complete ==="
    """
}


// -----------------------------------------------------------------------------
//  annotateAllDominantAlleles
// -----------------------------------------------------------------------------
/*process annotateAllDominantAlleles { 
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/all_dominant_annotations/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/all_dominant_annotations/${method}"

    input:
    tuple path(dominant_summary_tsv), val(method)

    output:
    path "${params.name_prefix}_all_dominant_annotated.tsv",         emit: annotated
    path "${params.name_prefix}_all_dominant_with_annotations.tsv",  emit: merged

    script:
    """
    python $projectDir/scripts/annotate_core_genes_sequential.py \\
        --input-summary ${dominant_summary_tsv} \\
        --output-annotations ${params.name_prefix}_all_dominant_annotated.tsv \\
        --output-merged ${params.name_prefix}_all_dominant_with_annotations.tsv \\
        --max-alternatives 30
    """
}*/


// -----------------------------------------------------------------------------
//  categorizeGOTerms
// -----------------------------------------------------------------------------
process categorizeGOTerms { 
    conda "$projectDir/envs/pangenome_env.yml"

    publishDir "${params.enhanced_output}/go_categorization/${method}", mode: 'copy'
    storeDir "${params.baseDir}/cache/go_categorization/${method}"

    input:
    tuple path(all_annotations), path(gene_matrix), path(gene_labels),
          path(core_genes), path(revigo_all_bp), path(revigo_all_mf),
          path(proteome_metadata), val(method)

    output:
    path "go_*_clusters.tsv"
    path "go_cluster_categorization_summary.txt"
    path "strain_cluster_analysis.tsv"
    path "*.png" optional true

    script:
    """
    echo "Categorizing GO terms with Revigo clusters..."

    python $projectDir/scripts/categorize_go_terms.py \\
        --all-annotations ${all_annotations} \\
        --gene-matrix ${gene_matrix} \\
        --gene-labels ${gene_labels} \\
        --core-genes ${core_genes} \\
        --revigo-all-bp ${revigo_all_bp} \\
        --revigo-all-mf ${revigo_all_mf} \\
        --proteome-metadata ${proteome_metadata} \\
        --output-dir .
    """
}
