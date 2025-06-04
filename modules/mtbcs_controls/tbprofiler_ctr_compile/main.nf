process tbprofiler_ctr_compile {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0.0
    @description: 
        This module runs TB-Profiler cimpile on the single sample using the TBDB database. 
        It first creates a symbolic link to the data.
    @changelog:
        v1.0.0-2024-11-01: Initial version
*/

    conda params.tbprofiler_env

        
    publishDir "${params.outDir}/${params.runID}/tbprofiler/", mode: 'copy'

    input:
        path(tbprofile_handover)


    output:
        path("tbprofiler.txt"),         emit: tbprofile_compiled
        path("tbprofiler.variants.csv")
        path("tbprofiler.variants.txt")
        

    script:
    """
    # Create sybolic links to the tbprofiler results
        ln -s ${params.outDir}/negative-controls/tbprofiler/results/ .
        ln -s ${params.outDir}/negative-controls/tbprofiler/bam/ .
        ln -s ${params.outDir}/negative-controls/tbprofiler/vcf/ .

    # Example command to compile input files
        tb-profiler collate
    """
}