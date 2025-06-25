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
        --samplesheet           [CSV]   Path to input data, must have the following columns (including the header):
                                            sampleID,forward_path,reverse_path,species,index
        --runID                 [str]   Run identifier
        --outDir                [path]  The output directory where the results will be saved
        --workDir               [path]  The temporary work directory for intermediate files
        --workflow              [str]   Workflow to run: full, comp_and_contam, tuberculosis_ctrls

    Optional arguments:
        --help                          Print this help message.

    Additional parameters:

    """
    exit 0
}

workflow {

    def color_green  = '\u001B[32m'
    def color_purp   = '\u001B[35m'
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

    // Show help if requested
    if (params.help) { 
        helpMessage() 
    }

    // Parameter validation
    def missingParams = []
    if (params.samplesheet == null) missingParams << "samplesheet"
    if (params.runID == null) missingParams << "runID"
    if (params.outDir == null) missingParams << "outDir"
    if (params.workDir == null) missingParams << "workDir"
    if (params.workflow == null) missingParams << "workflow"

    if (missingParams.size() > 0) {
        error "The following required parameters are missing: ${missingParams.join(', ')}. Please provide them with the appropriate flags."
        helpMessage()
    }

    // Validate workflow parameter
    def valid_workflows = ['full', 'comp_and_contam', 'tuberculosis_ctrls']
    if (!valid_workflows.contains(params.workflow)) {
        error "Invalid workflow parameter: ${params.workflow}. Valid options are: ${valid_workflows.join(', ')}"
        helpMessage()
    }

    /*
    ······································································································
        Main workflow - check taxonomic composition of reads and identify potential contamination
            - From the controls_ch, the samples are taxonomically classified with Sylph
            - Taxonomically classified sample reads and produces a summary of the reads
    ······································································································
    */

    if (params.workflow == 'full') {
        log.info "${color_cyan}Running the full workflow (${color_purp}--workflow full${color_cyan}) for taxonomic composition and contamination analysis.${color_reset}"
        comp_and_contam_wf(params.samplesheet)
        tuberculosis_wf(params.samplesheet)
        
    } else if (params.workflow == 'comp_and_contam') {
        log.info "${color_cyan}Running the single sample workflow (${color_purp}--workflow comp_and_contam${color_cyan}) for taxonomic composition and contamination analysis.${color_reset}"
        comp_and_contam_wf(params.samplesheet)
        
    } else if (params.workflow == 'tuberculosis_ctrls') {
        log.info "${color_cyan}Running the tuberculosis controls workflow (${color_purp}--workflow tuberculosis_ctrls${color_cyan}) for taxonomic composition and contamination analysis.${color_reset}"
        tuberculosis_wf(params.samplesheet)
    }

}