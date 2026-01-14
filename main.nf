#!/usr/bin/env nextflow

// Include processes
include { blast_db } from './modules/blast_db.nf'
include { blast_db_p } from './modules/blast_db_p.nf'
include { blast } from './modules/blast.nf'
include { blastp } from './modules/blastp.nf'
include { compare_lengths } from './modules/compare_lengths.nf'
include { concat_and_reformat } from './modules/concat_and_reformat.nf'
include { make_linking_file } from './modules/make_linking_file.nf'
include { reformat_bed } from './modules/reformat_bed.nf'
include { bedtools } from './modules/bedtools.nf'
include { transeq } from './modules/transeq.nf'
include { reformat_blast } from './modules/reformat_blast.nf'
include { identify_mut_prot } from './modules/identify_mut_prot.nf'
include { cd_hit } from './modules/cd_hit.nf'
include { identify_mut_clusters } from './modules/identify_mut_clusters.nf'
include { clustalo } from './modules/clustalo.nf'
include { identify_AAS } from './modules/identify_AAS.nf'
include { final_table } from './modules/final_table.nf'
include { final_mutants } from './modules/final_mutants.nf'
include { final_mutants_no_patterns } from './modules/final_mutants_no_patterns.nf'
include { graphing_r } from './modules/graphing_r.nf'
include { cluster_ani } from './modules/cluster_ani.nf'
include { ani_cluster_function } from './modules/ani_cluster_function.nf'

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
        if (params.run_main) {
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
            blast_db (
                host_genomes_ch
            )

            // Run blast of all input genes against all input genomes
            all_v_all = input_genes_ch
                .combine(blast_db.out.blast_db_path)

            blast(
                all_v_all
            )

            // Identify complete, incomplete, and missing genes
            fasta_blast_ch = blast.out
                .combine(input_genes_ch, by: 0)

            compare_lengths(
                fasta_blast_ch
            )

            // Obtain only gene ids
            all_gene_ids = input_genes_ch
                .map { gene_id, file -> gene_id }

            // Group the incomplete results
            incomplete_list_concat = concatGeneLists(all_gene_ids, (compare_lengths.out.Incomplete_list ?: Channel.empty()), "Incomplete")
            
            // Group the missing results
            missing_list_concat = concatGeneLists(all_gene_ids, (compare_lengths.out.Missing_list ?: Channel.empty()), "Missing")

            // Group the split results
            split_list_concat = concatGeneLists(all_gene_ids, (compare_lengths.out.Split_list ?: Channel.empty()), "Split")

            // Make linking file for downstream analysis
            make_linking_file(
                host_genomes_ch
            )
            
            collected_linking_files_ch =  make_linking_file.out
                .collectFile(
                    storeDir: "${params.outdir}/linking_file",
                    name: 'linking_file.tsv'
                )

            // Concatenate and reformat host genomes
            collected_fasta_files_ch = host_genomes_ch
                .map { it[1] }
                .collect()

            concat_and_reformat(
                collected_fasta_files_ch
            )

            // Concatenate complete blast results
            combined_blast = compare_lengths.out.Complete
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
            reformat_bed(
                combined_blast
            )

            // Run bedtools from blast result
            reformat_bed_fasta_ch = reformat_bed.out
                .combine( concat_and_reformat.out )

            bedtools(
                reformat_bed_fasta_ch
            )

            // Translate nucleotide to protein sequence
            transeq(
                bedtools.out
            )

            // Make a blast database of protein reference files
            blast_db_p(
                input_prot_ch
            )

            // Run blastp of protein reference files against translated sequences
            translated_prot_ch = transeq.out
                .combine(blast_db_p.out, by: 0)

            blastp(
                translated_prot_ch
            )

            // Filter blast results to identify WT and mutant sequences
            reformat_blast(
                blastp.out
            )

            // Filter for mutant protein sequences
            muts_accessions_prot_ch = reformat_blast.out.muts_accessions
                .combine( transeq.out, by: 0)

            identify_mut_prot(
                muts_accessions_prot_ch
            )

            // Cluster mutant protein sequences with cd-hit
            cd_hit(
                identify_mut_prot.out
            )

            // Filter for cd-hit reference mutants
            cluster_muts_list_fasta_ch = cd_hit.out.muts_list
                .combine( cd_hit.out.cd_hit_fasta, by:0 )

            identify_mut_clusters(
                cluster_muts_list_fasta_ch
            )

            split_fasta = identify_mut_clusters.out
                .splitFasta( file: true )
                .combine( input_prot_ch, by:0 )

            // Add reference proteins from PAO1 to all cd-hit reference sequences and align with clustalo
            clustalo(
                split_fasta
            )

            // Identify amino acid subtitution mutations
            aln_identify_AAS = clustalo.out
                .map { gene_id, file ->
                    def file_id = file.baseName.replace('_aln', '')
                    tuple(gene_id, file, file_id)
                }

            identify_AAS(
                aln_identify_AAS
            )

            // Combine alignment csvs

            combine_aln = identify_AAS.out
                .groupTuple()
                .map { gene_ID, files ->
                    def out_dir = file("${params.outdir}/mutant_calling_output/${gene_ID}/alignments")

                    def out_file = file("${out_dir}/${gene_ID}_combined_aln.csv")

                    out_file.text = files.collect { it.text }.join()

                    tuple(gene_ID, out_file)
                }
            
            // Generate final table
            combine_aln_split = combine_aln
                .combine( cd_hit.out.cd_hit_table, by:0 )

            final_table(
                combine_aln_split
            )

            // If true, generate presence absence files
            if (params.pres_abs) {

                if (!params.input_muts || params.input_muts == 'null') {
                    // Generate presence absence tables
                    incomplete_missing_ch = incomplete_list_concat
                        .combine( missing_list_concat, by:0 )
                        .combine( split_list_concat, by:0 )
                        .combine( final_table.out.final_mut_table, by:0 )
                        .combine( collected_linking_files_ch )
                        .combine( reformat_blast.out.muts_graphing, by:0 )

                    final_mutants_no_patterns(
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
                        .combine( final_table.out.final_mut_table, by:0 )
                        .combine( collected_linking_files_ch )
                        .combine( reformat_blast.out.muts_graphing, by:0 )
                        .combine( input_muts_ch, by:0 )

                    final_mutants(
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
                cluster_ani(
                    ani_values
                )

                collected_ani_clusters_ch =  cluster_ani.out
                    .collectFile(
                        storeDir: "${params.outdir}/combined_ani",
                        name: 'combined_ani_clusters.csv',
                        keepHeader: true
                    )

                if (params.input_muts == null) {
                    ani_functions_ch = final_mutants_no_patterns.out.functions_incomplete
                        .combine( collected_ani_clusters_ch )

                } else {
                    ani_functions_ch = final_mutants.out.functions_incomplete
                        .combine( collected_ani_clusters_ch )
                }

                ani_cluster_function(
                    ani_functions_ch
                )

            }

            // If true, graph phylogenetic trees and bargraphs based on gene presence absence files
            if (params.graphing) {

                graphing_input = final_mutants.out.functions_incomplete
                    .combine( final_mutants.out.functions_pres_abs_incomplete, by:0 )
                    .combine( ani_cluster_function.out.no_funct_clusters_env, by:0 )
                    
                graphing_r(
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
            
            no_funct_clusters_env = Channel
                .fromPath( params.no_funct_clusters )
                .map { file_path ->
                    def gene = file_path.parent.name
                    tuple(gene, file_path)
                }
                .groupTuple()

            graphing_input = functions_incomplete_ch
                .combine( functions_pres_abs_incomplete_ch, by:0 )
                .combine( no_funct_clusters_env, by:0 )
            
            graphing_r(
                graphing_input
            )
        
        }
    
}