#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.name_prefix            = 'Strep'

// ============================================================================
// DATABASE SELECTION: Choose annotation method
// ============================================================================
// For NCBI data: set params.database = "ncbi" and ensure GTF files are available
// For UniProt data: set params.database = "uniprot" (uses API, no GTF files needed)
// ============================================================================

// file paths
params.baseDir                = "$PWD"
//params.raw_input_directory    = "${params.baseDir}/simple_1258"
params.raw_input_directory    = "${params.baseDir}/fasta_foldseek"
params.gtf_directory          = "${params.baseDir}/gtf_1260"
//params.consolidated_dir       = "${params.baseDir}/1258_consol"
params.consolidated_dir       = "${params.baseDir}/foldseek_consolidated"
params.consolidated_fasta     = "${params.consolidated_dir}/consolidated.faa"
params.shared_headers_file    = "${params.consolidated_dir}/shared_headers.tsv"
params.metadata_file          = "${params.consolidated_dir}/metadata.tsv"
//params.proteome_metadata_file = "${params.baseDir}/streptomyces_genomes_metadata.tsv"
params.proteome_metadata_file = "${params.baseDir}/my_proteome_metadata.tsv"

// CD-HIT Revigo files
params.revigo_cdhit_bp_file   = "${params.baseDir}/revigo/cdhit_acc/Revigo_BP_Table.tsv"
params.revigo_cdhit_mf_file   = "${params.baseDir}/revigo/cdhit_acc/Revigo_MF_Table.tsv"

// Foldseek Revigo files
params.revigo_foldseek_bp_file = "${params.baseDir}/revigo/foldseek/Revigo_BP_Table.tsv"
params.revigo_foldseek_mf_file = "${params.baseDir}/revigo/foldseek/Revigo_MF_Table.tsv"

// SwiftOrtho Revigo files (for future use)
params.revigo_swiftortho_bp_file = "${params.baseDir}/revigo/swiftortho/Revigo_BP_Table.tsv"
params.revigo_swiftortho_mf_file = "${params.baseDir}/revigo/swiftortho/Revigo_MF_Table.tsv"

// parameters for all GO terms Revigo files
params.revigo_all_cdhit_bp_file = "${params.baseDir}/revigo/all_cdhit/Revigo_BP_Table.tsv"
params.revigo_all_cdhit_mf_file = "${params.baseDir}/revigo/all_cdhit/Revigo_MF_Table.tsv"

// Output directories
params.output_dir             = "${params.baseDir}/output"
params.enhanced_output        = "${params.output_dir}/uniprot_output"
params.cdhit_directory        = "${params.enhanced_output}/cdhit"
params.swiftortho_directory   = "${params.enhanced_output}/swiftortho"
params.foldseek_directory     = "${params.enhanced_output}/foldseek"
params.renamed_directory      = "${params.enhanced_output}/renamed_sequences"
params.sco_mapping            = "${params.baseDir}/SCO_NP_mapping.csv"

params.pdb_directory          = "${params.baseDir}/output/pdb_files/alphafold_pdbs_ncbi"
params.pdb_mapping            ="${params.baseDir}/output/pdb_files/final_pdb_map.tsv"

// Tools paths
params.swiftortho_path        = "/home/saba/SwiftOrtho"
params.foldseek_path          = '/usr/local/bin/foldseek'

// Clustering method selection
params.clustering_method      = "cdhit"  // Options: "cdhit", "swiftortho", "foldseek", "all"

// Database type selection for annotation
params.database               = "uniprot"  // Options: "ncbi", "uniprot"

// UniProt annotation optimization parameters
params.uniprot_batch_size     = 500        // Batch size for UniProt API requests (100-500)
params.uniprot_max_workers    = 8          // Concurrent workers for UniProt API (4-12)

// CD-HIT parameters
params.cdhit_identity         = 0.65
params.cdhit_coverage         = 0.75
params.threads                = 12
params.cdhit_mode            = "accurate"  // Options: "fast", "accurate"

// SwiftOrtho parameters
params.swiftortho_evalue      = "1e-5"
params.swiftortho_coverage    = 0.75
params.min_similarity         = 0.65

// Foldseek parameters
params.foldseek_sensitivity      = null    // Foldseek default: 9.5
params.foldseek_num_iterations   = null    // default: off
params.foldseek_max_seqs         = null    // default: 1000
params.foldseek_evalue           = null    // default: 0.001
params.foldseek_alignment_type   = null    // default: 2
params.foldseek_coverage         = null    // default: 0.0
params.foldseek_cov_mode         = null    // default: 0
params.foldseek_tmscore_threshold= null    // (no default)
params.foldseek_lddt_threshold   = null    // (no default)

// Parameters for using existing foldseek results
params.use_existing_foldseek_db      = true
params.existing_foldseek_db_path     = "${params.baseDir}/foldseek_done/db"
params.use_existing_foldseek_clusters = true
params.existing_foldseek_clusters_path = "${params.baseDir}/foldseek_done/clusters/clusters.tsv"

// how many CPUs to give to Foldseek (omit or set to null to use the default of 1)
params.foldseek_threads     = 10

// override only the minimum sequence identity (leave all other foldseek_* params null
// so they fall back to Foldseek's defaults)
params.foldseek_min_seq_id  = null
//parameter of missing core genes plot for scalibility
params.max_strains_missing_plot = 50  // Default: show top 50 strains in plot, 0 = show all

params.incremental_mode = false  // Enable incremental clustering
params.new_sequences_dir = null  // New sequences to add
params.previous_cdhit_repr = null  // Previous nr90 (representatives)
params.previous_cdhit_clstr = null  // Previous nr90.clstr

params.has_taxid = false  // Set to true if input headers are taxid|protein_id format
params.output_taxid = false  // Set to true for taxid|accession output, false for accession only


// Sequence database and annotation control
params.sequence_database = "uniprot"  // Options: "ncbi", "uniprot", "auto"
params.annotate_scope = "core"      // Options: "core", "all"
params.max_alternatives_per_cluster = 3  // Limit alternatives for efficiency



// Include shared processes
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

// Add these to the existing includes for SwiftOrtho
include { extractDominantAllelesWithCore as extractDominantAllelesWithCoreSwift } from './modules.nf'
include { annotateSmartDominantAlleles as annotateSmartDominantAllelesSwift } from './modules.nf'

include { analyzeGeneStructure as analyzeGeneStructureCDHIT }     from './modules.nf'
include { analyzeGeneStructure as analyzeGeneStructureSwift }    from './modules.nf'
include { analyzeGeneStructure as analyzeGeneStructureFoldseek } from './modules.nf'

// functional core analysis with Revigo

include { analyzeRevigoFunctionalCore as analyzeRevigoFunctionalCoreCDHIT }     from './modules.nf'
include { analyzeRevigoFunctionalCore as analyzeRevigoFunctionalCoreSwift }    from './modules.nf'
include { analyzeRevigoFunctionalCore as analyzeRevigoFunctionalCoreFoldseek } from './modules.nf'

include { clusterGOTermsWithGOATools as clusterGOTermsWithGOAToolsCDHIT } from './modules.nf'

// functional core analysis with GOATools clustering
include { analyzeClusteredFunctionalCore as analyzeClusteredFunctionalCoreCDHIT } from './modules.nf'

include { categorizeGOTerms as categorizeGOTermsCDHIT } from './modules.nf'

include { categorizeGOTerms as categorizeGOTermsFoldseek } from './modules.nf'

params.revigo_all_foldseek_bp_file = "${params.baseDir}/revigo/all_foldseek/Revigo_BP_Table.tsv"
params.revigo_all_foldseek_mf_file = "${params.baseDir}/revigo/all_foldseek/Revigo_MF_Table.tsv"

process consolidateSequences {
    publishDir "${params.consolidated_dir}", mode: 'copy'
    storeDir "${params.baseDir}/cache/consol"
    
    input:
    path input_dir
    
    output:
    path "consolidated.faa", emit: consolidated_fasta
    path "metadata.tsv", emit: metadata
    path "shared_headers.tsv", emit: shared_headers
    
    script:
    def has_taxid_arg = params.has_taxid ? 'true' : 'false'
    def output_taxid_arg = params.output_taxid ? 'true' : 'false'
    """
    python ${params.baseDir}/scripts/consolidate_sequences.py \\
        ${input_dir} \\
        consolidated.faa \\
        metadata.tsv \\
        shared_headers.tsv \\
        ${has_taxid_arg} \\
        ${output_taxid_arg}
    """
}

workflow {
    /*
     * Main entrypoint for the flexible pangenome pipeline
     * Works with any proteome data with taxid|protein_id headers
     */
    main:
        log.info "Starting pangenome pipeline"
        log.info "Base directory: ${params.baseDir}"
        log.info "Input directory: ${params.raw_input_directory}"

        // Check that required directories exist
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
        def doCDHIT    = methods.contains("cdhit")        || methods.contains("all")
        def doSwift    = methods.contains("swiftortho")   || methods.contains("all")
        def doFoldseek = methods.contains("foldseek")     || methods.contains("all")

        // Input channels
        input_dir_ch = Channel.fromPath(params.raw_input_directory)
        
        // Run consolidation if files don't exist
        if (!file(params.consolidated_fasta).exists() || 
            !file(params.shared_headers_file).exists() || 
            !file(params.metadata_file).exists()) {
            
            log.info "Consolidated files not found. Running consolidation..."
            consolidated_results = consolidateSequences(input_dir_ch)
            consolidated_fasta_ch = consolidated_results.consolidated_fasta
            shared_headers_ch = consolidated_results.shared_headers
        } else {
            log.info "Using existing consolidated files"
            consolidated_fasta_ch = Channel.fromPath(params.consolidated_fasta)
            shared_headers_ch = Channel.fromPath(params.shared_headers_file)
        }

        // Dummy mapping for non-foldseek renames
        def dummy = file("${params.enhanced_output}/dummy_mapping.tsv")
        dummy.text = "UniProt_ID\tPDB_file\tStatus\n"
        dummy_mapping_ch = Channel.fromPath(dummy.toString())

        if (doCDHIT) {
            def hasRevigoFiles = file(params.revigo_cdhit_bp_file).exists() && 
                                file(params.revigo_cdhit_mf_file).exists()
            
            if (!hasRevigoFiles) {
                log.warn "CD-HIT Revigo files not found - skipping Revigo-based functional analysis"
            }
            
            log.info "Running CD-HIT clustering method"
            cdhit_out = runCDHIT(consolidated_fasta_ch)
            clstr_ch = cdhit_out.map { it[0] }
            cdhit_input = clstr_ch.combine(consolidated_fasta_ch)

            renamed_cdhit = renameCDHITSequences(
                cdhit_input,         // tuple(cluster.clstr, fasta)
                'cdhit',             // method name
                shared_headers_ch,   // shared-headers.tsv
                dummy_mapping_ch     // dummy mapping
            )

            cdhit_tables = generateCDHITPangenomeTables(renamed_cdhit, input_dir_ch)
            validateCDHITPangenome(cdhit_tables, renamed_cdhit.map{ tuple(it[0], it[1]) }, input_dir_ch)
            
            // Core genome analysis
            core_genome_results = analyzeAndValidateCoreGenomeCDHIT(
                cdhit_tables.map{ tuple(it[1], it[3]) }, 
                'cdhit'
            )
            
            // Extract core genes channel from results
            core_genes_ch = core_genome_results.map { it[0] }  // Extract core_genes.txt
            
            // Combine for unified dominant allele extraction
            cdhit_combined = cdhit_tables.combine(renamed_cdhit.map{ tuple(it[0], it[1]) })
            cdhit_with_core = cdhit_combined.combine(core_genes_ch)
            
            // NEW: Unified dominant allele extraction with core separation
            unified_dominant_cdhit = extractDominantAllelesWithCoreCDHIT(cdhit_with_core)
            
            // Extract summary files for annotation
            annotation_input_ch = unified_dominant_cdhit.map { all_faa, all_summary, core_faa, core_summary, method ->
                tuple(all_summary, core_summary, method)
            }

            // Run annotation (new script: emits only *_core_genes_annotated.tsv and *_accessory_genes_annotated.tsv)
            annotation_results = annotateSmartDominantAllelesCDHIT(annotation_input_ch)
            core_annotations_ch = annotation_results.core_annotations
            accessory_annotations_ch = annotation_results.accessory_annotations
            
            // Run GO term clustering with GOATools
            /*
            goatools_input_ch = core_annotations_ch
                .combine(accessory_annotations_ch)
                .map { core_ann, accessory_ann -> tuple(core_ann, accessory_ann, 'cdhit') }
            
            goatools_results = clusterGOTermsWithGOAToolsCDHIT(goatools_input_ch)
*/
            // Create proteome metadata channel
            proteome_metadata_ch = Channel.fromPath(params.proteome_metadata_file, checkIfExists: true)

            // Extract gene matrix and labels from cdhit_tables
            gene_matrix_npz_ch = cdhit_tables.map { it[1] }   // strain_by_gene.npz
            gene_labels_ch     = cdhit_tables.map { it[3] }   // strain_by_gene.npz.labels.txt

            // Extract allele names for cluster mapping optimization
            allele_names_ch = renamed_cdhit.map { it[1] }  // Extract allele names TSV
            
            // GOATools-based Functional Core Analysis
            /*
            goatools_functional_input_cdhit_ch = gene_matrix_npz_ch
                .combine(gene_labels_ch)
                .combine(core_genes_ch)
                .combine(core_annotations_ch)
                .combine(proteome_metadata_ch)
                .combine(goatools_results.core_bp)
                .combine(goatools_results.core_mf)
                .combine(allele_names_ch)
                .map { matrix, labels, core_genes, core_ann, proteome_metadata, goatools_bp, goatools_mf, allele_names ->
                    tuple(matrix, labels, core_genes, core_ann, proteome_metadata, goatools_bp, goatools_mf, allele_names, 'cdhit')
                }*/

            //goatools_functional_results_cdhit = analyzeClusteredFunctionalCoreCDHIT(goatools_functional_input_cdhit_ch)

            // Heaps' law analysis (unchanged)
            heapslaw_cdhit_results = analyzeHeapsLawCDHIT(
                gene_matrix_npz_ch.combine(gene_labels_ch),
                "cdhit"
            )

            // Gene-Level Structural Analysis (now feed core annotations file)
            gene_structural_input_ch = gene_matrix_npz_ch
                .combine(gene_labels_ch)
                .combine(core_genes_ch)
                .combine(proteome_metadata_ch)
                .combine(core_annotations_ch)
                .map { matrix, labels, core_genes, proteome_metadata, core_ann ->
                    tuple(matrix, labels, core_genes, proteome_metadata, core_ann, 'cdhit')
                }

            structural_results_cdhit = analyzeGeneStructureCDHIT(gene_structural_input_ch)

            // Revigo analysis (if files exist) — feed core annotations file
 
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
                log.info "Running analyzeRevigoFunctionalCore anyway for testing - creating placeholder Revigo files"
                
                // Create placeholder Revigo files for testing
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
        
            //all_dominant_summary_ch = dominant_cdhit.map { faa, summary, method -> tuple(summary, method) }
            // all_dominant_annotations_ch = annotateAllDominantAllelesCDHIT(all_dominant_summary_ch)
            //all_dominant_results = annotateAllDominantAllelesCDHIT(all_dominant_summary_ch)

            //all_dominant_results = annotateAllDominantAllelesCDHIT(all_dominant_summary_ch)

            // Create channels for all GO terms Revigo files
            /*
            revigo_all_bp_ch = Channel.fromPath(params.revigo_all_cdhit_bp_file, checkIfExists: true)
            revigo_all_mf_ch = Channel.fromPath(params.revigo_all_cdhit_mf_file, checkIfExists: true)

            // NEW: Categorize GO terms with Revigo files
            go_categorization_input = all_dominant_results.merged
                .combine(gene_matrix_npz_ch)
                .combine(gene_labels_ch)
                .combine(core_genes_ch)
                .combine(revigo_all_bp_ch)
                .combine(revigo_all_mf_ch)
                .combine(proteome_metadata_ch)  // Add this line
                .map { annot, matrix, labels, core, bp, mf, metadata -> 
                    tuple(annot, matrix, labels, core, bp, mf, metadata, 'cdhit') 
                }

            go_categorization_results = categorizeGOTermsCDHIT(go_categorization_input)*/
        }
 
        // SwiftOrtho branch (uses consolidated FASTA)
        if (doSwift) {
            def hasRevigoFilesSwift = file(params.revigo_swiftortho_bp_file).exists() && 
                                    file(params.revigo_swiftortho_mf_file).exists()
            
            if (!hasRevigoFilesSwift) {
                log.warn "SwiftOrtho Revigo files not found - skipping Revigo-based functional analysis"
            }
            
            log.info "Running SwiftOrtho clustering method"
            swift_out = runSwiftOrtho(consolidated_fasta_ch)
            swift_clusters = swift_out.map{ it[2] }
            swift_input = swift_clusters.combine(consolidated_fasta_ch)

            renamed_swift = renameSwiftOrthoSequences(
                swift_input,
                'swiftortho',        // method name
                shared_headers_ch,
                dummy_mapping_ch
            )

            swift_tables = generateSwiftOrthoPangenomeTables(renamed_swift, input_dir_ch)
            validateSwiftOrthoPangenome(swift_tables, renamed_swift.map{ tuple(it[0], it[1]) }, input_dir_ch)
            
            // Core genome analysis
            core_genome_results_swift = analyzeAndValidateCoreGenomeSwift(
                swift_tables.map{ tuple(it[1], it[3]) }, 
                'swiftortho'
            )
            
            // Extract core genes channel from results
            core_genes_swift_ch = core_genome_results_swift.map { it[0] }  // Extract core_genes.txt
            
            // Combine for unified dominant allele extraction
            swift_combined = swift_tables.combine(renamed_swift.map{ tuple(it[0], it[1]) })
            swift_with_core = swift_combined.combine(core_genes_swift_ch)
            
            // Unified dominant allele extraction with core separation
            unified_dominant_swift = extractDominantAllelesWithCoreSwift(swift_with_core)
            
            // Extract summary files for annotation
            annotation_input_swift_ch = unified_dominant_swift.map { all_faa, all_summary, core_faa, core_summary, method ->
                tuple(all_summary, core_summary, method)
            }

            // Run annotation
            annotation_results_swift = annotateSmartDominantAllelesSwift(annotation_input_swift_ch)
            core_annotations_swift_ch = annotation_results_swift.core_annotations
            accessory_annotations_swift_ch = annotation_results_swift.accessory_annotations
            
            // Create proteome metadata channel
            proteome_metadata_ch = Channel.fromPath(params.proteome_metadata_file, checkIfExists: true)

            // Extract gene matrix and labels from swift_tables
            gene_matrix_npz_swift_ch = swift_tables.map { it[1] }   // strain_by_gene.npz
            gene_labels_swift_ch     = swift_tables.map { it[3] }   // strain_by_gene.npz.labels.txt

            // Extract allele names for cluster mapping optimization
            allele_names_swift_ch = renamed_swift.map { it[1] }  // Extract allele names TSV
            
            // Heaps' law analysis
            heapslaw_swift_results = analyzeHeapsLawSwiftOrtho(
                gene_matrix_npz_swift_ch.combine(gene_labels_swift_ch),
                "swiftortho"
            )

            // Gene-Level Structural Analysis (feed core annotations file)
            gene_structural_input_swift_ch = gene_matrix_npz_swift_ch
                .combine(gene_labels_swift_ch)
                .combine(core_genes_swift_ch)
                .combine(proteome_metadata_ch)
                .combine(core_annotations_swift_ch)
                .map { matrix, labels, core_genes, proteome_metadata, core_ann ->
                    tuple(matrix, labels, core_genes, proteome_metadata, core_ann, 'swiftortho')
                }

            structural_results_swift = analyzeGeneStructureSwift(gene_structural_input_swift_ch)

            // Revigo analysis (if files exist) — feed core annotations file
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
                // Optionally, you could create placeholder files like in CD-HIT branch for testing
            }
        }

        if (doFoldseek) {
            def hasRevigoFilesFoldseek = file(params.revigo_foldseek_bp_file).exists() && 
                                         file(params.revigo_foldseek_mf_file).exists()
            if (!hasRevigoFilesFoldseek) {
                log.warn "Foldseek Revigo files not found - skipping Revigo-based functional analysis"
            }

            log.info "Running Foldseek clustering method"

            // 1) Get clusters: use existing or build DB + cluster
            if (params.use_existing_foldseek_clusters && file(params.existing_foldseek_clusters_path).exists()) {
                log.info "Using existing Foldseek clusters from: ${params.existing_foldseek_clusters_path}"
                foldseek_clusters_ch = Channel.fromPath(params.existing_foldseek_clusters_path)
            } else {
                // Verify PDB directory
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

            // 2) Rename sequences with real PDB mapping
            foldseek_input = foldseek_clusters_ch.combine(consolidated_fasta_ch)
            pdb_mapping_ch = Channel.fromPath(params.pdb_mapping, checkIfExists: true)

            renamed_foldseek = renameFoldseekSequences(
                foldseek_input,
                'foldseek',
                shared_headers_ch,
                pdb_mapping_ch
            )

            // 3) Build pangenome tables and validate
            foldseek_tables = generateFoldseekPangenomeTables(renamed_foldseek, input_dir_ch)
            validateFoldseekPangenome(
                foldseek_tables,
                renamed_foldseek.map { tuple(it[0], it[1]) },
                input_dir_ch
            )

            // 4) Core-genome analysis
            core_genome_results_foldseek = analyzeAndValidateCoreGenomeFoldseek(
                foldseek_tables.map { tuple(it[1], it[3]) },  // (gene_matrix, gene_labels)
                'foldseek'
            )
            core_genes_foldseek_ch = core_genome_results_foldseek.map { it[0] } // core_genes.txt

            // 5) Dominant alleles (unified) + annotation
            // Combine: (allele_npz, gene_npz, allele_labels, gene_labels, method) ⨂ (renamed_faa, allele_names)
            foldseek_combined     = foldseek_tables.combine(renamed_foldseek.map { tuple(it[0], it[1]) })
            foldseek_with_core    = foldseek_combined.combine(core_genes_foldseek_ch)

            unified_dominant_foldseek = extractDominantAllelesWithCoreFoldseek(foldseek_with_core)

            // Annotation scope respected by process param
            annotation_input_foldseek_ch = unified_dominant_foldseek.map { all_faa, all_summary, core_faa, core_summary, method ->
                tuple(all_summary, core_summary, method)
            }
            annotation_results_foldseek   = annotateSmartDominantAllelesFoldseek(annotation_input_foldseek_ch)
            core_annotations_foldseek_ch  = annotation_results_foldseek.core_annotations
            accessory_annotations_foldseek_ch = annotation_results_foldseek.accessory_annotations

            // 6) Heaps’ law + structural analysis
            proteome_metadata_ch            = Channel.fromPath(params.proteome_metadata_file, checkIfExists: true)
            gene_matrix_npz_foldseek_ch     = foldseek_tables.map { it[1] } // strain_by_gene.npz
            gene_labels_foldseek_ch         = foldseek_tables.map { it[3] } // labels

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

            // 7) Revigo functional analysis (optional, only if files exist)
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

            
            //log.info "Foldseek comprehensive analysis completed"

            // Annotate ALL dominant alleles for Foldseek
            /*
            all_dominant_summary_foldseek_ch = dominant_foldseek.map { faa, summary, method -> tuple(summary, method) }
            all_dominant_results_foldseek = annotateAllDominantAllelesFoldseek(all_dominant_summary_foldseek_ch)

            // Create channels for all GO terms Revigo files (Foldseek)
            revigo_all_foldseek_bp_ch = Channel.fromPath(params.revigo_all_foldseek_bp_file, checkIfExists: true)
            revigo_all_foldseek_mf_ch = Channel.fromPath(params.revigo_all_foldseek_mf_file, checkIfExists: true)
                
            // Categorize GO terms with Revigo files for Foldseek
            go_categorization_input_foldseek = all_dominant_results_foldseek.merged
                .combine(gene_matrix_npz_foldseek_ch)
                .combine(gene_labels_foldseek_ch)
                .combine(core_genes_foldseek_ch)
                .combine(revigo_all_foldseek_bp_ch)
                .combine(revigo_all_foldseek_mf_ch)
                .combine(proteome_metadata_ch)
                .map { annot, matrix, labels, core, bp, mf, metadata -> 
                    tuple(annot, matrix, labels, core, bp, mf, metadata, 'foldseek') 
                }

            go_categorization_results_foldseek = categorizeGOTermsFoldseek(go_categorization_input_foldseek)*/
            
            //log.info "Foldseek comprehensive analysis with GO categorization completed"
}

process runCDHIT {
    publishDir "${params.cdhit_directory}", mode: 'copy'
    storeDir "${params.baseDir}/cache/cdhit"

    input:
    path non_redundant_seqs

    output:
    tuple path("cdhit_output.faa.clstr"), path("cdhit_output.faa")

    script:
    // Set CD-HIT mode based on parameter
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

process runSwiftOrtho {
    publishDir "${params.swiftortho_directory}/clusters", mode: 'copy'
    storeDir "${params.baseDir}/cache/swiftortho"
    
    input:
    path nr_faa_swift
    
    output:
    tuple path("input.fsa.sc"), path("input.fsa.sc.orth"), path("input.fsa.sc.orth.apc")

    script:
    """
    echo "SwiftOrtho path: ${params.swiftortho_path}"
    source /home/saba/anaconda3/etc/profile.d/conda.sh
    conda activate ${params.baseDir}/envs/swiftortho_env

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

process foldseekCreateDB {
  publishDir "${params.foldseek_directory}/database", mode: 'copy'
  storeDir   "${params.baseDir}/cache/foldseek_db"

  input:
    path pdb_dir

  output:
    path "foldseek_db"  // Output the complete database directory

  script:
  """
  # Create a new directory for the database
  mkdir -p foldseek_db
  
  # Create the database
  ${params.foldseek_path} createdb ${pdb_dir} foldseek_db/targetDB
  
  # Create the index (important for faster searches)
  ${params.foldseek_path} createindex foldseek_db/targetDB tmp_index
  
  # Verify database was created successfully
  echo "Database files:"
  ls -la foldseek_db/
  """
}

process foldseekCluster {
  publishDir "${params.foldseek_directory}/clusters", mode: 'copy'
  
  input:
    path db_dir
    
  output:
    path "clusters.tsv"

  script:
  // Only include parameters that are explicitly set (not null)
  def sensitivity_param = params.foldseek_sensitivity ? "-s ${params.foldseek_sensitivity}" : ""
  def num_iterations_param = params.foldseek_num_iterations ? "--num-iterations ${params.foldseek_num_iterations}" : ""
  def max_seqs_param = params.foldseek_max_seqs ? "--max-seqs ${params.foldseek_max_seqs}" : ""
  def evalue_param = params.foldseek_evalue ? "-e ${params.foldseek_evalue}" : ""
  def alignment_type_param = params.foldseek_alignment_type ? "--alignment-type ${params.foldseek_alignment_type}" : ""
  def coverage_param = params.foldseek_coverage ? "-c ${params.foldseek_coverage}" : ""
  def cov_mode_param = params.foldseek_cov_mode ? "--cov-mode ${params.foldseek_cov_mode}" : ""
  def tmscore_threshold_param = params.foldseek_tmscore_threshold ? "--tmscore-threshold ${params.foldseek_tmscore_threshold}" : ""
  def lddt_threshold_param = params.foldseek_lddt_threshold ? "--lddt-threshold ${params.foldseek_lddt_threshold}" : ""
  def threads_param = params.foldseek_threads ? "--threads ${params.foldseek_threads}" : ""
  def min_seq_id_param = params.foldseek_min_seq_id ? "--min-seq-id ${params.foldseek_min_seq_id}" : ""
  
  """
  # Check database directory contents
  echo "Database directory contents:"
  ls -la ${db_dir}
  
  # Run clustering - FIXED SYNTAX: foldseek cluster <inputDB> <outputDB> <tmpDir>
  ${params.foldseek_path} cluster ${db_dir}/targetDB cluster_result tmp \
    ${min_seq_id_param} \
    ${sensitivity_param} \
    ${num_iterations_param} \
    ${max_seqs_param} \
    ${evalue_param} \
    ${alignment_type_param} \
    ${coverage_param} \
    ${cov_mode_param} \
    ${tmscore_threshold_param} \
    ${lddt_threshold_param} \
    ${threads_param}
    
  # Convert the clustering results to TSV format - FIXED SYNTAX
  ${params.foldseek_path} createtsv ${db_dir}/targetDB cluster_result cluster_result.tsv
  
  # Check if the output file was created and rename it
  if [ -f cluster_result.tsv ]; then
    mv cluster_result.tsv clusters.tsv
  else
    echo "Clustering failed, no output file created" >&2
    echo "Contents of current directory:"
    ls -la
    exit 1
  fi
  """
}
