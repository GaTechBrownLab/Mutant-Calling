process make_linking_file {
    input:
    tuple val(gene_ID), path(fasta_file)

    output:
    path("${gene_ID}_linking.tsv")

    script:
    """
    grep '^>' ${fasta_file} | sed 's/^>//; s/ .*//' | awk -v id=${gene_ID} '{print id "\\t" \$0}' > ${gene_ID}_linking.tsv
    """
}