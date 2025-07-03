process reformat_bed {

    conda './envs/python_scripts.yaml'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
    tuple val(gene_ID), path(blast_table), path(scripts)

    output:
    tuple val(gene_ID), path("${gene_ID}/${gene_ID}_bed.bed"), emit: bed_file

    script:
    """
    mkdir -p ${gene_ID}
    python "${scripts}/reformat_bed.py" "$blast_table" "$gene_ID"
    """
}