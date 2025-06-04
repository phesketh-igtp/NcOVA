process merge_results {

    input:
    tuple val(sample_id), path(results)

    output:
    tuple val(sample_id), path("merged_results/${sample_id}.tsv")

    script:
    """
    mkdir -p merged_results
    cat ${results} > merged_results/${sample_id}.tsv
    """
}