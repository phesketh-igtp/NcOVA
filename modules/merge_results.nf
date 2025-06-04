process merge_results {

    conda params.general_env

    publishDir "${params.outDir}/${params.runID}/results/", mode: 'copy'

    input:
        path(ska_distance)
        path(sylph_taxonomy)

    output:
        path("${params.runID}_merged_plate_results.html"), emit: ska_merged_res
        path("${params.runID}_merged_sylph_results.tsv"), emit: sylph_merged_res

    script:
    """
    

    exit 1
    """
}