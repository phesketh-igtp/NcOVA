process merge_results {

    conda params.r_env

    publishDir "${params.outDir}/${params.runID}/", mode: 'copy'

    input:
        path(ska_distance, stageAs: "ska_distance.csv")
        path(sylph_taxonomy, stageAs: "sylph.csv")

    output:
        path("${params.runID}_merged_plate_results.html"), emit: ska_merged_res
        path("${params.runID}_merged_sylph_results.tsv"), emit: sylph_merged_res

    script:
    """
    Rscript ${params.r_script_dir}/merge_results.R

    exit 1
    """
}