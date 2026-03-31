process cluster_ani {
    label 'light'
    publishDir "${params.outdir}/cluster_ani_output", mode: 'copy'

    input:
    val level

    output:
    path "${level}_output.csv"

    script:
    """
    clusterANI.py -i "$params.fastani" -o "${level}_output.csv" -t "$level" -p ${task.cpus}
    
    """
}