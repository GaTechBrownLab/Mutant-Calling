#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mutant-Calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/GaTechBrownLab/Mutant-Calling
    Author: Iris Irby <iirby3@gatech.edu> <irisirby2018@gmail.com>
----------------------------------------------------------------------------------------
*/

//
// HELP MESSAGE
//
def helpMessage() {
    log.info"""
    =========================================
     Mutant-Calling
    =========================================
    Usage:
    nextflow run main.nf -with-conda --input_genes ./data/genes/*.fna --input_prots ./data/prot/*.faa
        --input_muts ./data/mutation_patterns/*.txt --host_genomes ./data/hosts/*.fna
        --outdir /results --pres_abs true --graphing false

    Required arguments:
        --help                      Print this help message and exit
        --input_genes               Path to reference genes of interest in .fna format (ex: PAO1 lasR)
        --input_prots               Path to reference proteins of interest in .faa format (ex: PAO1 lasR)
        --input_muts                Path to text file with mutational pattern of interest (ex: W60R)
                                    If there is no file, no mutational patterns will be taken into account
        --host_genomes              Path to bacterial genomes of interest to compare to the reference
                                    genes of interest.
        --outdir                    Output directory for results
        --run_main                  Run main mutant-calling script
                                    (default: true)
        --run_graphing              Re-run phylogenetic graphing after already previously running run_main
                                    (default: false)

    Blast thresholds and options:
        --tolerance                 When calling incomplete genes, allow for a tolerance of 50 base pairs
                                    to ensure insertions are called correctly
                                    (default: 50)
        --gene_proportion           Length cutoff to define a complete gene
                                    (gene length > gene length * (gene length - gene_proportion))
                                    (default: 0.05)
        --split_tolerance           Length in base pairs difference from contig length to determine split
                                    (default: 5)
        --gene_difference           Length cutoff to define a missing gene
                                    (gene length < gene length * gene_difference)
                                    (default: 0.4)
        --WT_cutoff                 Percentage cutoff for a gene variant to be considered wild-type
                                    (default: 100)
        --mut_cutoff                Percentage cutoff for a gene variant to be considered the correct gene
                                    Less than cutoff would indicate that the gene is incorrect, very mutated,
                                    or deleted
                                    (default: 40)

    Phylogenetic analysis and graphing parameters (INCOMPLETE):
        --pres_abs                  Generate presence absence files defining each genomes
                                    gene version as functional, non-functional, or potential alternate function
                                    (default: true)
        --graphing                  Graph the presence absence results on a phylogenetic tree
                                    (default: false)
        --tree                      Input phylogenetic tree if graphing is enabled
        --tree_root                 Input root for the phylogenetic tree if graphing is enabled
        --fastani                   Input FastANI of individuals in the tree if graphing is enabled
        --env_data                  Input associated environmental metadata if graphing is enabled
        --lineage_probability       The probability of a node's state being present must be over 0.9
                                    to be in a lineage in ancestral state reconstruction
    
    ANI clustering parameters:
        --ani_clusters              Run ANI clustering to define genome relatedness
                                    (default: false)
        --ani_start                 Start level for ANI clustering
                                    (default: 99.5)
        --ani_end                   End level for ANI clustering
                                    (default: 99.99)
        --ani_steps                 Steps between ani_start and ani_end
                                    (default: 15)
    
    run_graphing parameters:
        --incomplete_ran            Path for functions for graphing
        --pa_incomplete_ran         Path for presence absence matrix for graphing
        --no_funct_clusters         Path for non-functional clusters per environment

    Performance options:
        --blast_threads             Specify number of threads for blast
                                    (default: 10)
    Profiles:
        standard                    local execution
    
    Author:
        Iris Irby (iirby3@gatech.edu)
    """.stripIndent()
}


// Include processes
include { BLAST_DB } from './modules/BLAST_DB.nf'
include { BLAST_DB_P } from './modules/BLAST_DB_P.nf'
include { BLAST } from './modules/BLAST.nf'
include { BLASTP } from './modules/BLASTP.nf'
include { COMPARE_LENGTHS } from './modules/COMPARE_LENGTHS.nf'
include { CONCAT_AND_REFORMAT } from './modules/CONCAT_AND_REFORMAT.nf'
include { MAKE_LINKING_FILE } from './modules/MAKE_LINKING_FILE.nf'
include { REFORMAT_BED } from './modules/REFORMAT_BED.nf'
include { BEDTOOLS } from './modules/BEDTOOLS.nf'
include { TRANSEQ } from './modules/TRANSEQ.nf'
include { REFORMAT_BLAST } from './modules/REFORMAT_BLAST.nf'
include { IDENTIFY_MUT_PROT } from './modules/IDENTIFY_MUT_PROT.nf'
include { CDHIT } from './modules/CDHIT.nf'
include { IDENTIFY_MUT_CLUSTERS } from './modules/IDENTIFY_MUT_CLUSTERS.nf'
include { CLUSTALO } from './modules/CLUSTALO.nf'
include { IDENTIFY_AAS } from './modules/IDENTIFY_AAS.nf'
include { FINAL_TABLE } from './modules/FINAL_TABLE.nf'
include { FINAL_MUTANTS } from './modules/FINAL_MUTANTS.nf'
include { FINAL_MUTANTS_NO_PATTERNS } from './modules/FINAL_MUTANTS_NO_PATTERNS.nf'
include { GRAPHING_R } from './modules/GRAPHING_R.nf'
include { CLUSTER_ANI } from './modules/CLUSTER_ANI.nf'
include { ANI_CLUSTER_FUNCTION } from './modules/ANI_CLUSTER_FUNCTION.nf'

// Function for concatenating gene lists
def concatGeneLists(all_gene_ids, channel, suffix) {
    def grouped = channel
        .ifEmpty {
            def gene_ID = ""
            def files = ""

            tuple(gene_ID, files)
        }
        .groupTuple()

    return all_gene_ids
        .join(grouped, remainder: true)
            .map { gene_ID, files ->
                def out_dir = file("${params.outdir}/mutant_calling_output/${gene_ID}")
                // out_dir.mkdirs

                def out_file = file("${out_dir}/${suffix}_output.txt")

                if (gene_ID instanceof String && gene_ID.isEmpty()) {

                } else if (files && files.size() > 0) {
                    out_file.text = files.collect { it.text }.join()

                    tuple(gene_ID, out_file)
                } else {
                    out_file.text = ""

                    tuple(gene_ID, out_file)
                }
            }
}

// Begin main workflow
workflow {
    main:
        // Show help message
        if (params.help){
            helpMessage()
            exit 0
        }

        if (!params.run_main & !params.run_graphing) {
                exit 1, "No input provided! Please set `--run_main` or `run_graphing` to true. See `--help` for more details."
        }

        if (params.run_main) {

            // Validate parameters
            if (params.input_genes == "${projectDir}/data/genes/*.fna") {
                log.warn "The example data is being run, no `--input_genes` provided."
            }
            
            if (params.input_prots == "${projectDir}/data/prot/*.faa") {
                log.warn "The example data is being run, no `--input_prots` provided."
            }

            if (params.host_genomes == "${projectDir}/data/hosts/*.fna") {
                log.warn "The example data is being run, no `--host_genomes` provided."
            }

            if (params.input_muts == "${projectDir}/data/mutation_patterns/*.txt" ) {
                log.warn "The example data is being run, no `--input_muts` provided."
            }

            if (!params.input_muts) {
                log.warn "No mutants provided, the pipeline will be run without specific mutants."
            }

            // Set input channels
            input_genes_ch = Channel.fromPath( params.input_genes )
                    .map { file -> 
                    def id = file.baseName.replace('_nucl', '')
                    tuple(id, file)
                }
            
            input_prot_ch = Channel.fromPath( params.input_prots )
                    .map { file -> 
                    def id = file.baseName.replace('_prot', '')
                    tuple(id, file)
                }

            host_genomes_ch = Channel.fromPath( params.host_genomes )
                    .map { file -> tuple(file.baseName, file) }

            // Make blast db of all host genomes
            BLAST_DB (
                host_genomes_ch
            )

            // Run blast of all input genes against all input genomes
            all_v_all = input_genes_ch
                .combine(BLAST_DB.out.blast_db_path)

            BLAST(
                all_v_all
            )

            // Identify complete, incomplete, and missing genes
            fasta_blast_ch = BLAST.out
                .combine(input_genes_ch, by: 0)

            COMPARE_LENGTHS(
                fasta_blast_ch
            )

            // Obtain only gene ids
            all_gene_ids = input_genes_ch
                .map { gene_id, file -> gene_id }

            // Group the incomplete results
            incomplete_list_concat = concatGeneLists(all_gene_ids, (COMPARE_LENGTHS.out.Incomplete_list ?: Channel.empty()), "Incomplete")
            
            // Group the missing results
            missing_list_concat = concatGeneLists(all_gene_ids, (COMPARE_LENGTHS.out.Missing_list ?: Channel.empty()), "Missing")

            // Group the split results
            split_list_concat = concatGeneLists(all_gene_ids, (COMPARE_LENGTHS.out.Split_list ?: Channel.empty()), "Split")

            // Make linking file for downstream analysis
            MAKE_LINKING_FILE(
                host_genomes_ch
            )
            
            collected_linking_files_ch =  MAKE_LINKING_FILE.out
                .collectFile(
                    storeDir: "${params.outdir}/linking_file",
                    name: 'linking_file.tsv'
                )

            // Concatenate and reformat host genomes
            collected_fasta_files_ch = host_genomes_ch
                .map { it[1] }
                .collect()

            CONCAT_AND_REFORMAT(
                collected_fasta_files_ch
            )

            // Concatenate complete blast results
            combined_blast = COMPARE_LENGTHS.out.Complete
                .groupTuple()
                .map { gene_ID, files ->
                    def out_dir = file("${params.outdir}/mutant_calling_output/${gene_ID}")

                    def out_file = file("${out_dir}/${gene_ID}_complete_combined_blast.txt")

                    def combined_text = files.collect { f ->
                        def base = f.baseName.replaceAll("_${gene_ID}\$", "")
                        
                        f.readLines().collect { line ->
                            "${base}\t${line}"
                        }.join('\n')
                    }.join('\n')

                    out_file.text = combined_text

                    tuple(gene_ID, out_file)
                }
            
            // Create bed file from blast results
            REFORMAT_BED(
                combined_blast
            )

            // Run bedtools from blast result
            reformat_bed_fasta_ch = REFORMAT_BED.out
                .combine( CONCAT_AND_REFORMAT.out )

            BEDTOOLS(
                reformat_bed_fasta_ch
            )

            // Translate nucleotide to protein sequence
            TRANSEQ(
                BEDTOOLS.out
            )

            // Make a blast database of protein reference files
            BLAST_DB_P(
                input_prot_ch
            )

            // Run blastp of protein reference files against translated sequences
            translated_prot_ch = TRANSEQ.out
                .combine(BLAST_DB_P.out, by: 0)

            BLASTP(
                translated_prot_ch
            )

            // Filter blast results to identify WT and mutant sequences
            REFORMAT_BLAST(
                BLASTP.out
            )

            // Filter for mutant protein sequences
            muts_accessions_prot_ch = REFORMAT_BLAST.out.muts_accessions
                .combine( TRANSEQ.out, by: 0)

            IDENTIFY_MUT_PROT(
                muts_accessions_prot_ch
            )

            // Cluster mutant protein sequences with cd-hit
            CDHIT(
                IDENTIFY_MUT_PROT.out
            )

            // Filter for cd-hit reference mutants
            cluster_muts_list_fasta_ch = CDHIT.out.muts_list
                .combine( CDHIT.out.cd_hit_fasta, by:0 )

            IDENTIFY_MUT_CLUSTERS(
                cluster_muts_list_fasta_ch
            )

            split_fasta = IDENTIFY_MUT_CLUSTERS.out
                .splitFasta( file: true )
                .combine( input_prot_ch, by:0 )

            // Add reference proteins from PAO1 to all cd-hit reference sequences and align with clustalo
            CLUSTALO(
                split_fasta
            )

            // Identify amino acid subtitution mutations
            aln_identify_AAS = CLUSTALO.out
                .map { gene_id, file ->
                    def file_id = file.baseName.replace('_aln', '')
                    tuple(gene_id, file, file_id)
                }

            IDENTIFY_AAS(
                aln_identify_AAS
            )

            // Combine alignment csvs

            combine_aln = IDENTIFY_AAS.out
                .groupTuple()
                .map { gene_ID, files ->
                    def out_dir = file("${params.outdir}/mutant_calling_output/${gene_ID}/alignments")

                    def out_file = file("${out_dir}/${gene_ID}_combined_aln.csv")

                    out_file.text = files.collect { it.text }.join()

                    tuple(gene_ID, out_file)
                }
            
            // Generate final table
            combine_aln_split = combine_aln
                .combine( CDHIT.out.cd_hit_table, by:0 )

            FINAL_TABLE(
                combine_aln_split
            )

            // If true, generate presence absence files
            if (params.pres_abs) {

                if (!params.input_muts || params.input_muts == 'null') {
                    // Generate presence absence tables
                    incomplete_missing_ch = incomplete_list_concat
                        .combine( missing_list_concat, by:0 )
                        .combine( split_list_concat, by:0 )
                        .combine( FINAL_TABLE.out.final_mut_table, by:0 )
                        .combine( collected_linking_files_ch )
                        .combine( REFORMAT_BLAST.out.muts_graphing, by:0 )

                    FINAL_MUTANTS_NO_PATTERNS(
                        incomplete_missing_ch
                    )

                } else {
                    // Read in mutation patterns
                    input_muts_ch = Channel.fromPath( params.input_muts )
                        .map { file -> 
                        def id = file.baseName.replace('_mut_patterns', '')
                        tuple(id, file)
                        }

                    // Generate presence absence tables
                    incomplete_missing_ch = incomplete_list_concat
                        .combine( missing_list_concat, by:0 )
                        .combine( split_list_concat, by:0 )
                        .combine( FINAL_TABLE.out.final_mut_table, by:0 )
                        .combine( collected_linking_files_ch )
                        .combine( REFORMAT_BLAST.out.muts_graphing, by:0 )
                        .combine( input_muts_ch, by:0 )

                    FINAL_MUTANTS(
                        incomplete_missing_ch
                    )
                }
            }

            if (params.ani_clusters) {

                // ANI steps to channel
                ani_values = Channel.from((0..<params.ani_steps).collect { i ->
                        params.ani_start + i * (params.ani_end - params.ani_start) / (params.ani_steps - 1)
                    })
                    .map { String.format("%.8f", it) }

                // Cluster ANI values
                CLUSTER_ANI(
                    ani_values
                )

                collected_ani_clusters_ch =  CLUSTER_ANI.out
                    .collectFile(
                        storeDir: "${params.outdir}/combined_ani",
                        name: 'combined_ani_clusters.csv',
                        keepHeader: true
                    )

                if (params.input_muts == null) {
                    ani_functions_ch = FINAL_MUTANTS_NO_PATTERNS.out.functions_incomplete
                        .combine( collected_ani_clusters_ch )

                } else {
                    ani_functions_ch = FINAL_MUTANTS.out.functions_incomplete
                        .combine( collected_ani_clusters_ch )
                }

                ANI_CLUSTER_FUNCTION(
                    ani_functions_ch
                )

            }

            // If true, graph phylogenetic trees and bargraphs based on gene presence absence files
            if (params.graphing) {

                graphing_input = FINAL_MUTANTS.out.functions_incomplete
                    .combine( FINAL_MUTANTS.out.functions_pres_abs_incomplete, by:0 )
                    .combine( ANI_CLUSTER_FUNCTION.out.no_funct_clusters_env, by:0 )
                    
                GRAPHING_R(
                    graphing_input
                )
            }
        }

        if (params.run_graphing) {

            functions_incomplete_ch = Channel
                .fromPath( params.incomplete_ran )
                .filter { !it.name.contains("_all_") }
                .map { file_path ->
                    def gene = file_path.parent.name
                    tuple(gene, file_path)
                }
                .groupTuple()
            
            functions_pres_abs_incomplete_ch = Channel
                .fromPath( params.pa_incomplete_ran )
                .map { file_path ->
                    def gene = file_path.parent.name
                    tuple(gene, file_path)
                }
                .groupTuple()
            
            // no_funct_clusters_env = Channel
            //     .fromPath( params.no_funct_clusters )
            //     .map { file_path ->
            //         def gene = file_path.parent.name
            //         tuple(gene, file_path)
            //     }
            //     .groupTuple()

            graphing_input = functions_incomplete_ch
                .combine( functions_pres_abs_incomplete_ch, by:0 )
                .view()
            
            GRAPHING_R(
                graphing_input
            )
        
        }
    
}