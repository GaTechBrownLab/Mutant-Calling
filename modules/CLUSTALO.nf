process clustalo {
    label 'light'
    tag "${gene_ID}"
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/alignments/clustalo", mode: 'copy'

    input:
        tuple val(gene_ID), path(input_file), path(ref)

    output:
        tuple val(gene_ID), path("*_aln.fasta")
        
    script:
    """
    # Clean up input file
    new_file=\$(head -1 "$input_file")
    new_file_clean=\${new_file//>}.faa
    mv "${input_file}" "\${new_file_clean}"

    # Add PAO1 reference as the first input
    cat "$ref" <(echo) "\${new_file_clean}" > "\${new_file_clean%.*}_PAO1.faa"
    sed -i 's/ //g' "\${new_file_clean%.*}_PAO1.faa"

    # Run clustalo on cleaned input file to align
    clustalo -i \${new_file_clean%.*}_PAO1.faa -o \${new_file_clean%.*}_aln.fasta

    """
}