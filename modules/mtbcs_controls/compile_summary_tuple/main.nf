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

    tag "${runID}"
        
    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
    }
        
    publishDir "${params.outDir}/negative-controls/tbprofiler/", mode: 'copy'

    input:
        val(runID)
        val(complete_sampleID)


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