process COMPARE_LENGTHS {

    tag "COMPARE_LENGTHS: ${gene_ID}-${gene_file}"

    input:
        tuple val(gene_ID), path(table), path(gene_file)

    output:
        path "**"
        tuple val(gene_ID), path("${gene_ID}/Blast_gene_status/Complete/**"), optional: true, emit: Complete
        tuple val(gene_ID), path("${gene_ID}/Incomplete_output.txt"), optional: true, emit: Incomplete_list
        tuple val(gene_ID), path("${gene_ID}/Missing_output.txt"), optional: true, emit: Missing_list
        tuple val(gene_ID), path("${gene_ID}/Split_output.txt"), optional: true, emit: Split_list

    script:
    """
    # Make work directories for blast outputs after comparing lengths to reference
    mkdir -p ${gene_ID}/Blast_gene_status/Complete \
        ${gene_ID}/Blast_gene_status/Incomplete \
        ${gene_ID}/Blast_gene_status/Missing \
        ${gene_ID}/Blast_gene_status/Split
    
    # Obtain reference gene length
    length=\$(grep -v "^>" "$gene_file" | tr -d '\n' | wc -c)

    # Run compare_lengths.py
    compare_lengths.py \
        "\$length" \
        "$table" \
        "${gene_ID}" \
        "$params.gene_proportion" \
        "$params.gene_difference" \
        "$params.tolerance" \
        "$params.split_tolerance"

"""
}