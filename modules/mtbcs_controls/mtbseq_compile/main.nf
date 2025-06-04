process mtbseq_ctr_compile {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-03
    @version: 0.1
    @description: 
        This module compiles on the single sample using the TBDB database. 
        The module loops through the negative control outputs and captures all 
        the MTBSeq results and creates a new tabular files with the results
    @changelog:
        v1.0.1-2025-04-08: Change - Changed from loop to straight concatenate
*/

    tag "${runID}"
        
    conda params.mtbseq_env
        
    publishDir "${params.outDir}/${params.runID}/mtbseq/", mode: 'copy'

    input:
        val(runID)
        val(complete_sampleID)

    output:
        path("Strain_Classification.tab"),          emit: mtbseq_class_compiled
        path("Mapping_and_Variant_Statistics.tab"), emit: mtbseq_stats_compiled

    script:
    """
    # Concatenate the files
        cat ${params.outDir}/negative-controls/mtbseq/Classification/*.Strain_Classification.tab \\
                | sed '/^Date/d' | sed "s@'@@g" \\
                > Strain_Classification.tab

        cat ${params.outDir}/negative-controls/mtbseq/Statistics/*.Mapping_and_Variant_Statistics.tab \\
                | sed '/^Date/d' | sed "s@'@@g" \\
                > Mapping_and_Variant_Statistics.tab
    """
}