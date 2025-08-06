process cd_hit {
    
    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(prot_fasta)

    output:
        tuple val(gene_ID), path("${gene_ID}/*.clstr"), optional: true, emit: cd_hit_cluster
        tuple val(gene_ID), path("${gene_ID}/${gene_ID}_muts_cd_100.faa"), optional:true, emit: cd_hit_fasta

    script:
    """
    mkdir -p ${gene_ID}

    cd-hit -i $prot_fasta -o ${gene_ID}/${gene_ID}_muts_cd_100.faa -c 1.00 -n 5 -s 1 -d 200
    
    """
}

