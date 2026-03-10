#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
//   1. params block REMOVED from here entirely.
//      All params now live in nextflow.config where they belong.
//
//   2. consolidateSequences process:
//      ${params.baseDir}/input_scripts/  →  $projectDir/input_scripts/
//      $projectDir always points to where main.nf lives, regardless of who runs it or from where.
//
//   3. runSwiftOrtho process:
//      SwiftOrtho needs its own env because it uses a different Python env.
// ============================================================================

// Include shared processes from modules.nf
include { renameSequences as renameCDHITSequences }        from './modules.nf'
include { renameSequences as renameSwiftOrthoSequences }   from './modules.nf'
include { renameSequences as renameFoldseekSequences }     from './modules.nf'
include { generatePangenomeTables as generateCDHITPangenomeTables }   from './modules.nf'
include { generatePangenomeTables as generateSwiftOrthoPangenomeTables } from './modules.nf'
include { generatePangenomeTables as generateFoldseekPangenomeTables }   from './modules.nf'
include { validatePangenome as validateCDHITPangenome }     from './modules.nf'
include { validatePangenome as validateSwiftOrthoPangenome } from './modules.nf'
include { validatePangenome as validateFoldseekPangenome }   from './modules.nf'
include { viewPangenome as viewCDHITPangenome }            from './modules.nf'
include { viewPangenome as viewSwiftOrthoPangenome }       from './modules.nf'
include { viewPangenome as viewFoldseekPangenome }         from './modules.nf'
include { analyzeAndValidateCoreGenome as analyzeAndValidateCoreGenomeCDHIT } from './modules.nf'
include { analyzeAndValidateCoreGenome as analyzeAndValidateCoreGenomeSwift } from './modules.nf'
include { analyzeAndValidateCoreGenome as analyzeAndValidateCoreGenomeFoldseek } from './modules.nf'

include { analyzeHeapsLaw as analyzeHeapsLawCDHIT } from './modules.nf'
include { analyzeHeapsLaw as analyzeHeapsLawSwiftOrtho } from './modules.nf'
include { analyzeHeapsLaw as analyzeHeapsLawFoldseek } from './modules.nf'

include { extractDominantAllelesWithCore as extractDominantAllelesWithCoreCDHIT } from './modules.nf'
include { annotateSmartDominantAlleles as annotateSmartDominantAllelesCDHIT } from './modules.nf'

include { extractDominantAllelesWithCore as extractDominantAllelesWithCoreFoldseek } from './modules.nf'
include { annotateSmartDominantAlleles as annotateSmartDominantAllelesFoldseek } from './modules.nf'

include { extractDominantAllelesWithCore as extractDominantAllelesWithCoreSwift } from './modules.nf'
include { annotateSmartDominantAlleles as annotateSmartDominantAllelesSwift } from './modules.nf'

include { analyzeGeneStructure as analyzeGeneStructureCDHIT }     from './modules.nf'
include { analyzeGeneStructure as analyzeGeneStructureSwift }    from './modules.nf'
include { analyzeGeneStructure as analyzeGeneStructureFoldseek } from './modules.nf'

include { analyzeRevigoFunctionalCore as analyzeRevigoFunctionalCoreCDHIT }     from './modules.nf'
include { analyzeRevigoFunctionalCore as analyzeRevigoFunctionalCoreSwift }    from './modules.nf'
include { analyzeRevigoFunctionalCore as analyzeRevigoFunctionalCoreFoldseek } from './modules.nf'

include { clusterGOTermsWithGOATools as clusterGOTermsWithGOAToolsCDHIT } from './modules.nf'

include { analyzeClusteredFunctionalCore as analyzeClusteredFunctionalCoreCDHIT } from './modules.nf'

include { categorizeGOTerms as categorizeGOTermsCDHIT } from './modules.nf'
include { categorizeGOTerms as categorizeGOTermsFoldseek } from './modules.nf'


// ============================================================================
// consolidateSequences
// ============================================================================
process consolidateSequences {
    publishDir "${params.consolidated_dir}", mode: 'copy'
    storeDir "${params.baseDir}/cache/consol"
    conda "$projectDir/envs/pangenome_env.yml"

    input:
    path input_dir

    output:
    path "consolidated.faa",    emit: consolidated_fasta
    path "metadata.tsv",        emit: metadata
    path "shared_headers.tsv",  emit: shared_headers

    script:
    def has_taxid_arg  = params.has_taxid  ? 'true' : 'false'
    def output_taxid_arg = params.output_taxid ? 'true' : 'false'
    """
    python $projectDir/input_scripts/consolidate_sequences.py \\
        ${input_dir} \\
        consolidated.faa \\
        metadata.tsv \\
        shared_headers.tsv \\
        ${has_taxid_arg} \\
        ${output_taxid_arg}
    """
}


// ============================================================================
// runCDHIT
// conda directive added so cd-hit is guaranteed to be available from the env.
// ============================================================================
process runCDHIT {
    publishDir "${params.cdhit_directory}", mode: 'copy'
    storeDir "${params.baseDir}/cache/cdhit"

    // cd-hit is declared in pangenome_env.yml (bioconda package)
    conda "$projectDir/envs/pangenome_env.yml"

    input:
    path non_redundant_seqs

    output:
    tuple path("cdhit_output.faa.clstr"), path("cdhit_output.faa")

    script:
    def cdhit_mode_flag = params.cdhit_mode == "accurate" ? "-g 1" : "-g 0"
    """
    cd-hit -i ${non_redundant_seqs} \\
        -o cdhit_output.faa \\
        -c ${params.cdhit_identity} \\
        -n 4 \\
        -aL ${params.cdhit_coverage} \\
        -G 1 \\
        ${cdhit_mode_flag} \\
        -T ${params.threads} \\
        -M 10000
    """
}


// ============================================================================
// runSwiftOrtho
//   SwiftOrtho has its own env because it needs a specific Python + BLAST setup
//   separate from pangenome_env.
// ============================================================================
process runSwiftOrtho {
    publishDir "${params.swiftortho_directory}/clusters", mode: 'copy'
    storeDir "${params.baseDir}/cache/swiftortho"

    // CHANGE: conda directive replaces the manual activation inside the script
    conda "$projectDir/envs/swiftortho_env.yml"

    input:
    path nr_faa_swift

    output:
    tuple path("input.fsa.sc"), path("input.fsa.sc.orth"), path("input.fsa.sc.orth.apc")

    script:
    // CHANGE: removed `source detect_conda.sh` and `conda activate/deactivate`
    // Everything else is identical to original
    """
    echo "SwiftOrtho path: ${params.swiftortho_path}"

    python ${params.swiftortho_path}/bin/find_hit.py \\
        -p blastp \\
        -i ${nr_faa_swift} \\
        -d ${nr_faa_swift} \\
        -o input.fsa.sc \\
        -e ${params.swiftortho_evalue} \\
        -s 110101101111

    python ${params.swiftortho_path}/bin/find_orth.py \\
        -i input.fsa.sc \\
        -c ${params.min_similarity} \\
        -y ${params.swiftortho_coverage} \\
        -s '|' > input.fsa.sc.orth

    python ${params.swiftortho_path}/bin/find_cluster.py \\
        -i input.fsa.sc.orth \\
        -a apc \\
        -I 1 > input.fsa.sc.orth.apc
    """
}


// ============================================================================
// foldseekCreateDB
// ============================================================================
process foldseekCreateDB {
    publishDir "${params.foldseek_directory}/database", mode: 'copy'
    storeDir   "${params.baseDir}/cache/foldseek_db"

    input:
    path pdb_dir

    output:
    path "foldseek_db"

    script:
    """
    mkdir -p foldseek_db

    ${params.foldseek_path} createdb ${pdb_dir} foldseek_db/targetDB

    ${params.foldseek_path} createindex foldseek_db/targetDB tmp_index

    echo "Database files:"
    ls -la foldseek_db/
    """
}


// ============================================================================
// foldseekCluster
// ============================================================================
process foldseekCluster {
    publishDir "${params.foldseek_directory}/clusters", mode: 'copy'

    input:
    path db_dir

    output:
    path "clusters.tsv"

    script:
    def sensitivity_param      = params.foldseek_sensitivity      ? "-s ${params.foldseek_sensitivity}"                  : ""
    def num_iterations_param   = params.foldseek_num_iterations   ? "--num-iterations ${params.foldseek_num_iterations}"  : ""
    def max_seqs_param         = params.foldseek_max_seqs         ? "--max-seqs ${params.foldseek_max_seqs}"              : ""
    def evalue_param           = params.foldseek_evalue           ? "-e ${params.foldseek_evalue}"                        : ""
    def alignment_type_param   = params.foldseek_alignment_type   ? "--alignment-type ${params.foldseek_alignment_type}"  : ""
    def coverage_param         = params.foldseek_coverage         ? "-c ${params.foldseek_coverage}"                      : ""
    def cov_mode_param         = params.foldseek_cov_mode         ? "--cov-mode ${params.foldseek_cov_mode}"              : ""
    def tmscore_threshold_param= params.foldseek_tmscore_threshold? "--tmscore-threshold ${params.foldseek_tmscore_threshold}" : ""
    def lddt_threshold_param   = params.foldseek_lddt_threshold   ? "--lddt-threshold ${params.foldseek_lddt_threshold}"  : ""
    def threads_param          = params.foldseek_threads          ? "--threads ${params.foldseek_threads}"                : ""
    def min_seq_id_param       = params.foldseek_min_seq_id       ? "--min-seq-id ${params.foldseek_min_seq_id}"          : ""

    """
    echo "Database directory contents:"
    ls -la ${db_dir}

    ${params.foldseek_path} cluster ${db_dir}/targetDB cluster_result tmp \\
        ${min_seq_id_param} \\
        ${sensitivity_param} \\
        ${num_iterations_param} \\
        ${max_seqs_param} \\
        ${evalue_param} \\
        ${alignment_type_param} \\
        ${coverage_param} \\
        ${cov_mode_param} \\
        ${tmscore_threshold_param} \\
        ${lddt_threshold_param} \\
        ${threads_param}

    ${params.foldseek_path} createtsv ${db_dir}/targetDB cluster_result cluster_result.tsv

    if [ -f cluster_result.tsv ]; then
        mv cluster_result.tsv clusters.tsv
    else
        echo "Clustering failed, no output file created" >&2
        ls -la
        exit 1
    fi
    """
}


// ============================================================================
// WORKFLOW
// ============================================================================
workflow {
    main:
        log.info "Starting pangenome pipeline"
        log.info "Base directory: ${params.baseDir}"
        log.info "Input directory: ${params.raw_input_directory}"

        if (!file(params.raw_input_directory).exists()) {
            error "Input directory does not exist: ${params.raw_input_directory}"
        }

        // Create output directories
        file(params.output_dir).mkdirs()
        file(params.enhanced_output).mkdirs()
        file(params.cdhit_directory).mkdirs()
        file(params.swiftortho_directory).mkdirs()
        file(params.foldseek_directory).mkdirs()
        file(params.renamed_directory).mkdirs()
        file(params.consolidated_dir).mkdirs()

        // Decide which methods to run
        def methods    = params.clustering_method.split(",")
        def doCDHIT    = methods.contains("cdhit")      || methods.contains("all")
        def doSwift    = methods.contains("swiftortho") || methods.contains("all")
        def doFoldseek = methods.contains("foldseek")   || methods.contains("all")

        // Input channels
        input_dir_ch = Channel.fromPath(params.raw_input_directory)

        // Run consolidation if files don't exist
        if (!file(params.consolidated_fasta).exists() ||
            !file(params.shared_headers_file).exists() ||
            !file(params.metadata_file).exists()) {

            log.info "Consolidated files not found. Running consolidation..."
            consolidated_results  = consolidateSequences(input_dir_ch)
            consolidated_fasta_ch = consolidated_results.consolidated_fasta
            shared_headers_ch     = consolidated_results.shared_headers
        } else {
            log.info "Using existing consolidated files"
            consolidated_fasta_ch = Channel.fromPath(params.consolidated_fasta)
            shared_headers_ch     = Channel.fromPath(params.shared_headers_file)
        }

        // Dummy mapping for non-foldseek renames
        def dummy = file("${params.enhanced_output}/dummy_mapping.tsv")
        dummy.text = "UniProt_ID\tPDB_file\tStatus\n"
        dummy_mapping_ch = Channel.fromPath(dummy.toString())


        // ====================================================================
        // CD-HIT branch
        // ====================================================================
        if (doCDHIT) {
            def hasRevigoFiles = file(params.revigo_cdhit_bp_file).exists() &&
                                 file(params.revigo_cdhit_mf_file).exists()

            if (!hasRevigoFiles) {
                log.warn "CD-HIT Revigo files not found - skipping Revigo-based functional analysis"
            }

            log.info "Running CD-HIT clustering method"
            cdhit_out   = runCDHIT(consolidated_fasta_ch)
            clstr_ch    = cdhit_out.map { it[0] }
            cdhit_input = clstr_ch.combine(consolidated_fasta_ch)

            renamed_cdhit = renameCDHITSequences(
                cdhit_input,
                'cdhit',
                shared_headers_ch,
                dummy_mapping_ch
            )

            cdhit_tables = generateCDHITPangenomeTables(renamed_cdhit, input_dir_ch)
            validateCDHITPangenome(cdhit_tables, renamed_cdhit.map{ tuple(it[0], it[1]) }, input_dir_ch)

            // Core genome analysis
            core_genome_results = analyzeAndValidateCoreGenomeCDHIT(
                cdhit_tables.map{ tuple(it[1], it[3]) },
                'cdhit'
            )
            core_genes_ch = core_genome_results.map { it[0] }

            // Dominant allele extraction
            cdhit_combined  = cdhit_tables.combine(renamed_cdhit.map{ tuple(it[0], it[1]) })
            cdhit_with_core = cdhit_combined.combine(core_genes_ch)

            unified_dominant_cdhit = extractDominantAllelesWithCoreCDHIT(cdhit_with_core)

            annotation_input_ch = unified_dominant_cdhit.map { all_faa, all_summary, core_faa, core_summary, method ->
                tuple(all_summary, core_summary, method)
            }

            annotation_results        = annotateSmartDominantAllelesCDHIT(annotation_input_ch)
            core_annotations_ch       = annotation_results.core_annotations
            accessory_annotations_ch  = annotation_results.accessory_annotations

            // Proteome metadata channel
            proteome_metadata_ch = Channel.fromPath(params.proteome_metadata_file, checkIfExists: true)

            // Gene matrix channels
            gene_matrix_npz_ch = cdhit_tables.map { it[1] }
            gene_labels_ch     = cdhit_tables.map { it[3] }
            allele_names_ch    = renamed_cdhit.map { it[1] }

            // Heaps' law analysis
            heapslaw_cdhit_results = analyzeHeapsLawCDHIT(
                gene_matrix_npz_ch.combine(gene_labels_ch),
                "cdhit"
            )

            // Gene-level structural analysis
            gene_structural_input_ch = gene_matrix_npz_ch
                .combine(gene_labels_ch)
                .combine(core_genes_ch)
                .combine(proteome_metadata_ch)
                .combine(core_annotations_ch)
                .map { matrix, labels, core_genes, proteome_metadata, core_ann ->
                    tuple(matrix, labels, core_genes, proteome_metadata, core_ann, 'cdhit')
                }

            structural_results_cdhit = analyzeGeneStructureCDHIT(gene_structural_input_ch)

            // Revigo analysis (if files exist)
            if (hasRevigoFiles) {
                revigo_cdhit_bp_ch = Channel.fromPath(params.revigo_cdhit_bp_file)
                revigo_cdhit_mf_ch = Channel.fromPath(params.revigo_cdhit_mf_file)

                revigo_functional_input_cdhit_ch = gene_matrix_npz_ch
                    .combine(gene_labels_ch)
                    .combine(core_genes_ch)
                    .combine(core_annotations_ch)
                    .combine(proteome_metadata_ch)
                    .combine(revigo_cdhit_bp_ch)
                    .combine(revigo_cdhit_mf_ch)
                    .map { matrix, labels, core_genes, core_ann, proteome_metadata, revigo_bp, revigo_mf ->
                        tuple(matrix, labels, core_genes, core_ann, proteome_metadata, revigo_bp, revigo_mf, 'cdhit')
                    }

                revigo_functional_results_cdhit = analyzeRevigoFunctionalCoreCDHIT(revigo_functional_input_cdhit_ch)
            } else {
                log.info "Running analyzeRevigoFunctionalCore with placeholder Revigo files for testing"

                revigo_cdhit_bp_ch = Channel.of("Name\tTermID\tFrequency\tValue\tUniqueness\tDispensability\tRepresentative\nplaceholder\tGO:0000001\t1\t1.0\t1.0\t0.0\t1")
                    .collectFile(name: 'placeholder_bp.tsv', newLine: true, storeDir: "${params.baseDir}/temp")
                revigo_cdhit_mf_ch = Channel.of("Name\tTermID\tFrequency\tValue\tUniqueness\tDispensability\tRepresentative\nplaceholder\tGO:0000001\t1\t1.0\t1.0\t0.0\t1")
                    .collectFile(name: 'placeholder_mf.tsv', newLine: true, storeDir: "${params.baseDir}/temp")

                revigo_functional_input_cdhit_ch = gene_matrix_npz_ch
                    .combine(gene_labels_ch)
                    .combine(core_genes_ch)
                    .combine(core_annotations_ch)
                    .combine(proteome_metadata_ch)
                    .combine(revigo_cdhit_bp_ch)
                    .combine(revigo_cdhit_mf_ch)
                    .map { matrix, labels, core_genes, core_ann, proteome_metadata, revigo_bp, revigo_mf ->
                        tuple(matrix, labels, core_genes, core_ann, proteome_metadata, revigo_bp, revigo_mf, 'cdhit')
                    }

                revigo_functional_results_cdhit = analyzeRevigoFunctionalCoreCDHIT(revigo_functional_input_cdhit_ch)
            }
        }


        // ====================================================================
        // SwiftOrtho branch
        // ====================================================================
        if (doSwift) {
            def hasRevigoFilesSwift = file(params.revigo_swiftortho_bp_file).exists() &&
                                      file(params.revigo_swiftortho_mf_file).exists()

            if (!hasRevigoFilesSwift) {
                log.warn "SwiftOrtho Revigo files not found - skipping Revigo-based functional analysis"
            }

            log.info "Running SwiftOrtho clustering method"
            swift_out      = runSwiftOrtho(consolidated_fasta_ch)
            swift_clusters = swift_out.map{ it[2] }
            swift_input    = swift_clusters.combine(consolidated_fasta_ch)

            renamed_swift = renameSwiftOrthoSequences(
                swift_input,
                'swiftortho',
                shared_headers_ch,
                dummy_mapping_ch
            )

            swift_tables = generateSwiftOrthoPangenomeTables(renamed_swift, input_dir_ch)
            validateSwiftOrthoPangenome(swift_tables, renamed_swift.map{ tuple(it[0], it[1]) }, input_dir_ch)

            core_genome_results_swift = analyzeAndValidateCoreGenomeSwift(
                swift_tables.map{ tuple(it[1], it[3]) },
                'swiftortho'
            )
            core_genes_swift_ch = core_genome_results_swift.map { it[0] }

            swift_combined  = swift_tables.combine(renamed_swift.map{ tuple(it[0], it[1]) })
            swift_with_core = swift_combined.combine(core_genes_swift_ch)

            unified_dominant_swift = extractDominantAllelesWithCoreSwift(swift_with_core)

            annotation_input_swift_ch = unified_dominant_swift.map { all_faa, all_summary, core_faa, core_summary, method ->
                tuple(all_summary, core_summary, method)
            }

            annotation_results_swift       = annotateSmartDominantAllelesSwift(annotation_input_swift_ch)
            core_annotations_swift_ch      = annotation_results_swift.core_annotations
            accessory_annotations_swift_ch = annotation_results_swift.accessory_annotations

            proteome_metadata_ch = Channel.fromPath(params.proteome_metadata_file, checkIfExists: true)

            gene_matrix_npz_swift_ch = swift_tables.map { it[1] }
            gene_labels_swift_ch     = swift_tables.map { it[3] }
            allele_names_swift_ch    = renamed_swift.map { it[1] }

            heapslaw_swift_results = analyzeHeapsLawSwiftOrtho(
                gene_matrix_npz_swift_ch.combine(gene_labels_swift_ch),
                "swiftortho"
            )

            gene_structural_input_swift_ch = gene_matrix_npz_swift_ch
                .combine(gene_labels_swift_ch)
                .combine(core_genes_swift_ch)
                .combine(proteome_metadata_ch)
                .combine(core_annotations_swift_ch)
                .map { matrix, labels, core_genes, proteome_metadata, core_ann ->
                    tuple(matrix, labels, core_genes, proteome_metadata, core_ann, 'swiftortho')
                }

            structural_results_swift = analyzeGeneStructureSwift(gene_structural_input_swift_ch)

            if (hasRevigoFilesSwift) {
                revigo_swift_bp_ch = Channel.fromPath(params.revigo_swiftortho_bp_file)
                revigo_swift_mf_ch = Channel.fromPath(params.revigo_swiftortho_mf_file)

                revigo_functional_input_swift_ch = gene_matrix_npz_swift_ch
                    .combine(gene_labels_swift_ch)
                    .combine(core_genes_swift_ch)
                    .combine(core_annotations_swift_ch)
                    .combine(proteome_metadata_ch)
                    .combine(revigo_swift_bp_ch)
                    .combine(revigo_swift_mf_ch)
                    .map { matrix, labels, core_genes, core_ann, proteome_metadata, revigo_bp, revigo_mf ->
                        tuple(matrix, labels, core_genes, core_ann, proteome_metadata, revigo_bp, revigo_mf, 'swiftortho')
                    }

                revigo_functional_results_swift = analyzeRevigoFunctionalCoreSwift(revigo_functional_input_swift_ch)
            } else {
                log.info "Skipping Revigo functional analysis for SwiftOrtho - Revigo files not found"
            }
        }


        // ====================================================================
        // Foldseek branch
        // ====================================================================
        if (doFoldseek) {
            def hasRevigoFilesFoldseek = file(params.revigo_foldseek_bp_file).exists() &&
                                         file(params.revigo_foldseek_mf_file).exists()
            if (!hasRevigoFilesFoldseek) {
                log.warn "Foldseek Revigo files not found - skipping Revigo-based functional analysis"
            }

            log.info "Running Foldseek clustering method"

            // Get clusters: use existing or build from scratch
            if (params.use_existing_foldseek_clusters && file(params.existing_foldseek_clusters_path).exists()) {
                log.info "Using existing Foldseek clusters from: ${params.existing_foldseek_clusters_path}"
                foldseek_clusters_ch = Channel.fromPath(params.existing_foldseek_clusters_path)
            } else {
                def pdbDir = file(params.pdb_directory)
                if (!pdbDir.exists() || !pdbDir.isDirectory()) {
                    error "PDB directory does not exist: ${params.pdb_directory}"
                }
                def pdbFiles = pdbDir.listFiles().findAll { it.name.endsWith('.pdb') }
                if (pdbFiles.size() == 0) {
                    error "No PDB files found in directory: ${params.pdb_directory}"
                } else {
                    log.info "Found ${pdbFiles.size()} PDB files in ${params.pdb_directory}"
                }

                pdb_dir_ch = Channel.fromPath(params.pdb_directory)

                if (params.use_existing_foldseek_db && file(params.existing_foldseek_db_path).exists()) {
                    log.info "Using existing Foldseek database from: ${params.existing_foldseek_db_path}"
                    db_dir_ch = Channel.fromPath(params.existing_foldseek_db_path)
                } else {
                    db_dir_ch = foldseekCreateDB(pdb_dir_ch)
                }

                foldseek_clusters_ch = foldseekCluster(db_dir_ch)
            }

            // Rename with real PDB mapping
            foldseek_input  = foldseek_clusters_ch.combine(consolidated_fasta_ch)
            pdb_mapping_ch  = Channel.fromPath(params.pdb_mapping, checkIfExists: true)

            renamed_foldseek = renameFoldseekSequences(
                foldseek_input,
                'foldseek',
                shared_headers_ch,
                pdb_mapping_ch
            )

            // Pangenome tables and validation
            foldseek_tables = generateFoldseekPangenomeTables(renamed_foldseek, input_dir_ch)
            validateFoldseekPangenome(
                foldseek_tables,
                renamed_foldseek.map { tuple(it[0], it[1]) },
                input_dir_ch
            )

            // Core-genome analysis
            core_genome_results_foldseek = analyzeAndValidateCoreGenomeFoldseek(
                foldseek_tables.map { tuple(it[1], it[3]) },
                'foldseek'
            )
            core_genes_foldseek_ch = core_genome_results_foldseek.map { it[0] }

            // Dominant alleles + annotation
            foldseek_combined  = foldseek_tables.combine(renamed_foldseek.map { tuple(it[0], it[1]) })
            foldseek_with_core = foldseek_combined.combine(core_genes_foldseek_ch)

            unified_dominant_foldseek = extractDominantAllelesWithCoreFoldseek(foldseek_with_core)

            annotation_input_foldseek_ch = unified_dominant_foldseek.map { all_faa, all_summary, core_faa, core_summary, method ->
                tuple(all_summary, core_summary, method)
            }
            annotation_results_foldseek       = annotateSmartDominantAllelesFoldseek(annotation_input_foldseek_ch)
            core_annotations_foldseek_ch      = annotation_results_foldseek.core_annotations
            accessory_annotations_foldseek_ch = annotation_results_foldseek.accessory_annotations

            // Heaps' law + structural analysis
            proteome_metadata_ch            = Channel.fromPath(params.proteome_metadata_file, checkIfExists: true)
            gene_matrix_npz_foldseek_ch     = foldseek_tables.map { it[1] }
            gene_labels_foldseek_ch         = foldseek_tables.map { it[3] }

            heapslaw_foldseek_results = analyzeHeapsLawFoldseek(
                gene_matrix_npz_foldseek_ch.combine(gene_labels_foldseek_ch),
                "foldseek"
            )

            gene_structural_input_foldseek_ch = gene_matrix_npz_foldseek_ch
                .combine(gene_labels_foldseek_ch)
                .combine(core_genes_foldseek_ch)
                .combine(proteome_metadata_ch)
                .combine(core_annotations_foldseek_ch)
                .map { matrix, labels, core_genes, proteome_metadata, annotations ->
                    tuple(matrix, labels, core_genes, proteome_metadata, annotations, 'foldseek')
                }

            structural_results_foldseek = analyzeGeneStructureFoldseek(gene_structural_input_foldseek_ch)

            // Revigo functional analysis (optional)
            if (hasRevigoFilesFoldseek) {
                revigo_foldseek_bp_ch = Channel.fromPath(params.revigo_foldseek_bp_file)
                revigo_foldseek_mf_ch = Channel.fromPath(params.revigo_foldseek_mf_file)

                revigo_functional_input_foldseek_ch = gene_matrix_npz_foldseek_ch
                    .combine(gene_labels_foldseek_ch)
                    .combine(core_genes_foldseek_ch)
                    .combine(core_annotations_foldseek_ch)
                    .combine(proteome_metadata_ch)
                    .combine(revigo_foldseek_bp_ch)
                    .combine(revigo_foldseek_mf_ch)
                    .map { matrix, labels, core_genes, annotations, proteome_metadata, revigo_bp, revigo_mf ->
                        tuple(matrix, labels, core_genes, annotations, proteome_metadata, revigo_bp, revigo_mf, 'foldseek')
                    }

                revigo_functional_results_foldseek = analyzeRevigoFunctionalCoreFoldseek(revigo_functional_input_foldseek_ch)
            }
        }
}
