process clustalo {

    conda './envs/mutant_calling_genomics.yaml'

    publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(input_file), path(ref)

    output:
        tuple val(gene_ID), path("${gene_ID}/alignments/clustalo/*_aln.fasta")
        
    script:
    """
    mkdir -p ${gene_ID}/alignments/clustalo

    new_file=\$(head -1 "$input_file")
    new_file_clean=\${new_file//>}.faa
    mv "${input_file}" "\${new_file_clean}"

    cat "$ref" <(echo) "\${new_file_clean}" > "\${new_file_clean%.*}_PAO1.faa"
    sed -i 's/ //g' "\${new_file_clean%.*}_PAO1.faa"

    clustalo -i \${new_file_clean%.*}_PAO1.faa -o ${gene_ID}/alignments/clustalo/\${new_file_clean%.*}_aln.fasta

    """
}