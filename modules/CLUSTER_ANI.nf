process CLUSTER_ANI {

    tag "CLUSTER_ANI: ${level}"
    
    publishDir "${params.outdir}/cluster_ani_output", mode: 'copy'

    input:
    val level

    output:
    path "${level}_output.csv"

    script:
    """
    # Run clusterANI.py to cluster inputs into genomovars
    clusterANI.py -i "$params.fastani" -o "${level}_output.csv" -t "$level" -p 24
    
    """
}