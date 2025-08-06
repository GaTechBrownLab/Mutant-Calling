#!/usr/bin/env nextflow

// Include processes
include { blast_db } from './modules/blast_db.nf'
include { blast_db_p } from './modules/blast_db_p.nf'
include { blast } from './modules/blast.nf'
include { blastp } from './modules/blastp.nf'
include { compare_lengths } from './modules/compare_lengths.nf'
include { concat_and_reformat } from './modules/concat_and_reformat.nf'
include { combine_blast } from './modules/combine_blast.nf'
include { make_linking_file } from './modules/make_linking_file.nf'
include { reformat_bed } from './modules/reformat_bed.nf'
include { bedtools } from './modules/bedtools.nf'
include { transeq } from './modules/transeq.nf'
include { reformat_blast } from './modules/reformat_blast.nf'
include { identify_mut_prot } from './modules/identify_mut_prot.nf'
include { cd_hit } from './modules/cd_hit.nf'
include { split_files } from './modules/split_files.nf'
include { identify_mut_clusters } from './modules/identify_mut_clusters.nf'
include { clustalo } from './modules/clustalo.nf'
include { identify_AAS } from './modules/identify_AAS.nf'
include { combine_aln } from './modules/combine_aln.nf'
include { final_table } from './modules/final_table.nf'
include { final_mutants } from './modules/final_mutants.nf'
include { graphing_r } from './modules/graphing_r.nf'

// Begin main workflow
workflow {
    main:
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

        // Concatenate blast results
        grouped_complete_ch =  compare_lengths.out.Complete
            .groupTuple()

        combine_blast{
            grouped_complete_ch
        }
        
        // Create bed file from blast results
        reformat_bed(
            combine_blast.out
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

        // Reformat cd-hit final table
        split_files(
            cd_hit.out.cd_hit_cluster
        )

        // Filter for cd-hit reference mutants
        cluster_muts_list_fasta_ch = split_files.out.muts_list
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
        grouped_complete_ch =  identify_AAS.out
            .groupTuple()

        combine_aln(
           grouped_complete_ch 
        )
        
        // Generate final table
        combine_aln_split = combine_aln.out
            .combine( split_files.out.cd_hit_table, by:0 )

        final_table(
            combine_aln_split
        )

        // If true, generate presence absence files
        if (params.pres_abs) {
            // Read in mutation patterns
            input_muts_ch = Channel.fromPath( params.input_muts )
                .map { file -> 
                def id = file.baseName.replace('_mut_patterns', '')
                tuple(id, file)
            }
            
            // Generate list of genomes with incomplete gene
            grouped_incomplete_ch = compare_lengths.out.Incomplete_list
                .groupTuple()

            // Generate list of genomes with missing gene
            grouped_missing_ch = compare_lengths.out.Missing_list
                .groupTuple()

            // Generate presence absence tables
            incomplete_missing_ch = grouped_incomplete_ch
                .combine( grouped_missing_ch, by:0 )
                .combine( final_table.out.final_mut_table, by:0 )
                .combine( collected_linking_files_ch )
                .combine( reformat_blast.out.muts_graphing, by:0 )
                .combine( input_muts_ch, by:0 )

            final_mutants(
                incomplete_missing_ch
            )
        }

        // If true, graph phylogenetic trees and bargraphs based on gene presence absence files
        if (params.graphing) {

            graphing_input = final_mutants.out.functions_incomplete
                .combine( final_mutants.out.functions_pres_abs_incomplete, by:0 )
            
            graphing_r(
                graphing_input
            )
        }

}