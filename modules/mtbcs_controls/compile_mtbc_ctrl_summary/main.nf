process compile_mtbc_ctrl_summary {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-03
    @version: 1.2.0
    @description:
        This module concatenates all the statistical files within the output directory
        capturing all the results to date and not just in this run
    @changelog:
        v1.0.0-2025-04-03:
            Added - initial version
        v1.1.0-2025-04-04:
            Added - taxonkit to get formatted taxonomy
            Changed - output from TSV to CSV
        v1.2.0-2025-04-10:
            Changed - from 'taxonkit lineage' to 'taxonkit reformat' for consistency
*/

    conda params.taxonkit_env

    publishDir "${params.outDir}/${params.runID}/", mode: 'copy'

    input:
        path(mtbseq_class_collected)
        path(mtbseq_stats_collected)
        path(tbprofile_compiled)

    script:
    """
    exit 1
    """
}