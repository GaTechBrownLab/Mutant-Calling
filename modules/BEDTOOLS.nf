process BEDTOOLS {
    tag "${gene_ID}"
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}", mode: 'copy'

    input:
        tuple val(gene_ID), path(bed_file), path(fasta)

    output:
        tuple val(gene_ID), path("${gene_ID}_nucleotide_seq_cleaned.fna"), emit: bed_output_nucl

    script:
    """  
    # Bedtools
    bedtools getfasta -s -fi $fasta -bed $bed_file -fo ${gene_ID}_nucleotide_seq.fna

    # Clean output
    sed '/^>/ s/:.*//' ${gene_ID}_nucleotide_seq.fna > ${gene_ID}_nucleotide_seq_cleaned.fna

    """
}

