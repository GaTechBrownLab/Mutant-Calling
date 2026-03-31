process identify_mut_prot {
    tag "${gene_ID}"
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/cd_hit_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(muts_accessions), path(prot_fasta)

    output:
        tuple val(gene_ID), path("${gene_ID}_muts.faa"), optional: true, emit: mut_prots

    script:
    """
    mkdir -p ${gene_ID}
    
    awk '
    NR==FNR {
        ids[\$0];
        next
    } 
    /^>/ {
        header = \$0;
        sub(/^>/, "", header);  # Remove the leading ">" for the ID
        f = (header in ids);
    } 
    f' "$muts_accessions" "$prot_fasta" >> "${gene_ID}_muts.faa"

    sed -i 's/\\*/X/g' "${gene_ID}_muts.faa"

"""
}