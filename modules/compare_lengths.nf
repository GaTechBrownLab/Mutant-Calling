process compare_lengths {

        publishDir "${params.outdir}/mutant_calling_output", mode: 'copy'

    input:
        tuple val(gene_ID), path(table), path(gene_file)
        val tolerance
        val gene_proportion
        val gene_difference

    output:
        tuple val(gene_ID), path("${gene_ID}/Complete/**"), optional: true, emit: Complete
        tuple val(gene_ID), path("${gene_ID}/Incomplete/**"), optional: true, emit: Incomplete
        tuple val(gene_ID), path("${gene_ID}/Missing/**"), optional: true, emit: Missing
        tuple val(gene_ID), path("${gene_ID}/Status/Complete/**"), emit: Complete_list
        tuple val(gene_ID), path("${gene_ID}/Status/Incomplete/**"), emit: Incomplete_list
        tuple val(gene_ID), path("${gene_ID}/Status/Missing/**"), emit: Missing_list

    script:
    """
    mkdir -p ${gene_ID}
    mkdir -p ${gene_ID}/Complete ${gene_ID}/Incomplete ${gene_ID}/Missing

    mkdir -p ${gene_ID}/Status
    mkdir -p ${gene_ID}/Status/Complete ${gene_ID}/Status/Incomplete ${gene_ID}/Status/Missing

    length=\$(sed -n '2p' ${gene_file} | tr -d '\n' | wc -c)

    echo "Length: \$length"
    # Loop through all TSV files in the current directory

    touch ${gene_ID}/Status/Complete/${table}_complete.txt
    touch ${gene_ID}/Status/Incomplete/${table}_incomplete.txt
    touch ${gene_ID}/Status/Missing/${table}_missing.txt
        
    if [ ! -s "$table" ]; then
        echo "$table"
        # File is empty, move to Missing
        mv "$table" ${gene_ID}/Missing/
        echo "$table" >> ${gene_ID}/Status/Missing/${table}_missing.txt
    else
        # Count the number of lines in the file
        line_count=\$(wc -l < "$table")

        if [ "\$line_count" -eq 1 ]; then
            # File has only one line, move to Complete
            mv "$table" ${gene_ID}/Complete/
            echo "$table" >> ${gene_ID}/Status/Complete/${table}_complete.txt
        else
            # Read the first line and check columns 5 and 6
            first_line=\$(head -n 1 "$table")
            col7=\$(echo "\$first_line" | awk -F'\t' '{print \$7}')
            col6=\$(echo "\first_line" | awk -F'\t' '{print \$6}')

            second_line=\$(sed -n '2p' "$table")
            second_col7=\$(echo "\$second_line" | awk -F'\t' '{print \$7}')

            row1_col10=\$(echo "\$first_line" | awk -F'\t' '{print \$10}')
            row1_col11=\$(echo "\$first_line" | awk -F'\t' '{print \$11}')
            row2_col10=\$(echo "\$second_line" | awk -F'\t' '{print \$10}')
            row2_col11=\$(echo "\$second_line" | awk -F'\t' '{print \$11}')

            if [[ "\$col7" == "\$col6" ]]; then
                mv "$table" ${gene_ID}/Complete/
                echo "$table" >> ${gene_ID}/Status/Complete/${table}_complete.txt
            else
                length_threshold=\$(echo "\$length - \$length * $gene_proportion" | bc | cut -d'.' -f1)

                if [ "\$col7" -ge "\$length_threshold" ] && [ "\$col6" = "\$length" ]; then
                    mv "$table" ${gene_ID}/Complete/
                    echo "$table" >> ${gene_ID}/Status/Complete/${table}_complete.txt
                else
                    gene_diff_threshold=\$(echo "\$length * $gene_difference" | bc | cut -d'.' -f1)

                    if [ "\$col7" -le \$gene_diff_threshold ]; then
                        mv "$table" ${gene_ID}/Missing/
                        echo "$table" >> ${gene_ID}/Status/Missing/${table}_missing.txt
                    else
                        if [[ "\$row1_col10" -eq 1 ]]; then
                            x="\$row1_col11"
                            col10_other="\$row2_col10"
                            col11_other="\$row2_col11"
 
                            # Compute absolute differences
                            diff_col10=\$(( col10_other - x ))
                            diff_col11=\$(( col11_other - length ))
                            
                            if (( diff_col10 >= -$tolerance && diff_col10 <= $tolerance )) && \
                                (( diff_col11 >= -$tolerance && diff_col11 <= $tolerance )); then
                                mv "$table" ${gene_ID}/Incomplete/
                                echo "$table" >> ${gene_ID}/Status/Incomplete/${table}_incomplete.txt
                            else
                                mv "$table" ${gene_ID}/Complete/
                                echo "$table" >> ${gene_ID}/Status/Complete/${table}_complete.txt
                            fi

                        elif [[ "\$row2_col10" -eq 1 ]]; then
                            x="\$row2_col11"
                            col10_other="\$row1_col10"
                            col11_other="\$row1_col11"

                            # Compute absolute differences
                            diff_col10=\$(( col10_other - x ))
                            diff_col11=\$(( col11_other - length ))

                            if (( diff_col10 >= -$tolerance && diff_col10 <= $tolerance )) && \
                                (( diff_col11 >= -$tolerance && diff_col11 <= $tolerance )); then
                                mv "$table" ${gene_ID}/Incomplete/
                                echo "$table" >> ${gene_ID}/Status/Incomplete/${table}_incomplete.txt
                            else
                                mv "$table" ${gene_ID}/Complete/
                                echo "$table" >> ${gene_ID}/Status/Complete/${table}_complete.txt
                            fi
                        else
                            mv "$table" ${gene_ID}/Complete/
                            echo "$table" >> ${gene_ID}/Status/Complete/${table}_complete.txt
                        fi
                    fi
                fi
            fi
        fi
    fi

"""
}