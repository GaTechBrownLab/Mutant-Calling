process CDHIT {

    tag "CDHIT: ${gene_ID}"
    
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/cd_hit_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(prot_fasta)

    output:
        tuple val(gene_ID), path("*.clstr"), optional: true, emit: cd_hit_cluster
        tuple val(gene_ID), path("${gene_ID}_muts_cd_100.faa"), optional:true, emit: cd_hit_fasta
        tuple val(gene_ID), path("${gene_ID}_cd_hit_table.csv"), emit: cd_hit_table
        tuple val(gene_ID), path("${gene_ID}_muts.txt"), emit: muts_list

    script:
    """
    # Run cd-hit
    cd-hit -i $prot_fasta -o ${gene_ID}_muts_cd_100.faa -c 1.00 -n 5 -s 1 -d 200

    # Run reformat_cd_hit.py to reformat cd-hit
    reformat_cd_hit.py "$gene_ID" "${gene_ID}_muts_cd_100.faa.clstr"

    """
}

