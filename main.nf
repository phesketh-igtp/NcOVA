#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { comp_and_contam_wf }   from './workflow/comp_and_contam.nf'
include { tuberculosis_wf   }    from './workflow/tuberculosis_ctrls.nf'


/* 
    Help Message
*/

def helpMessage() {
    log.info """
    Usage:

    Mandatory arguments:
        --samplesheet           [CSV]   Path to input data (must be surrounded with quotes)
        --outDir                [path]  The output directory where the results will be saved
        --workDir               [path]  The temporary work directory for intermediate files (can be deleted when 
                                            analysis is complete to recovered storage space)

    Optional arguments:
        --help                          Print this help message.

    Additional parameters:


    """
    exit 0
}

/* 
    MAIN WORKFLOW
*/

workflow {

    //def color_purple = '\u001B[35m'
    def color_green  = '\u001B[32m'
    def color_red    = '\u001B[31m'
    def color_reset  = '\u001B[0m'
    def color_cyan   = '\u001B[36m'

    log.info """
    ${color_cyan}
    ════════════════════════════════════════════════════════════════════════
        ▗▖  ▗▖      ▗▄▖ ▗▖  ▗▖ ▗▄▖ 
        ▐▛▚▖▐▌ ▗▄▄▖▐▌ ▐▌▐▌  ▐▌▐▌ ▐▌
        ▐▌ ▝▜▌▐▌   ▐▌ ▐▌▐▌  ▐▌▐▛▀▜▌
        ▐▌  ▐▌▝▚▄▄▖▝▚▄▞▘ ▝▚▞▘ ▐▌ ▐▌                
        ${color_green}Pre-release development version${color_cyan}   
    ════════════════════════════════════════════════════════════════════════
        Negative cOntrol Validation Analysis ${params.version}
    ════════════════════════════════════════════════════════════════════════
    ${color_reset}
    """

    if (params.help) { helpMessage() }

    def missingParams = []
        if (params.samplesheet == null) missingParams << "samplesheet"
        if (params.runID == null) missingParams << "runID"
        if (params.outDir == null) missingParams << "outDir"
        if (params.workDir == null) missingParams << "workDir"

        if (missingParams.size() > 0) {
            error "The following required parameters are missing: ${missingParams.join(', ')}. Please provide them with the appropriate flags."
            helpMessage()
        }

        /*
            DEFINE INPUT ARGUMENTS: expected argument to be provided at time of running 
                nextflow at CLI
            nextflow run main.nf \
                --samplesheet /path/to/sample-sheet
                --runID [a-zA-Z0-9]
                --workflow [full, single, pairwise, summary, barcoding]
        */

        if (params.samplesheet == null) { error "Please provide a samplesheet CSV file with --samplesheet (csv). The sample sheet but have the following headers - sampleID,forward_path,reverse_path,kingdom,index"; helpMessage() }
        if (params.runID == null) { error "Please provide a runID file with --runID (chr)"; helpMessage() }
        if (params.outDir == null) { error "Please provide a results/database directory for the RutiSeq db (location where new or past results will be) with --outDir (path)"; helpMessage() }
        if (params.workDir == null) { error "Please provide a work directory for the temporary intermediate files --workDir (path)"; helpMessage() }


        /*
        ······································································································
            Main workflow - check taxonomic composition of reads and identify potential contamination
                - From the controls_ch, the samples are taxonomically classified with Sylph
                - Taxonomically classified sample reads and produces a summary of the reads
        ······································································································
        */

        comp_and_contam_wf( params.samplesheet )



}

/*
author: Poppy J Hesketh Best
date: 2025-04-04
version: 1.0.0-beta
description: 
    This is the main workflow for the RutiSeq-nf pipeline. It is designed to be run with Nextflow and 
        takes a samplesheet as input. The workflow performs the following steps:
            - Update the TBProfiler database
            - Perform negative control analysis
            - Perform single sample analysis
            - Perform pairwise sample analysis
            - Produce summary tables and visualisations
            - Perform barcoding analysis (optional-WIP)
changelog
        - 2024-11-01: Initial version
*/


