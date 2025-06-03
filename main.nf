#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { sylph_read_taxonomy   } from './modules/sylph_read_taxonomy.nf'
include { seqkit_stats          } from './modules/seqkit_stats.nf'
include { ska_reads_build       } from './modules/read_ska.nf'
include { ska_distance          } from './modules/ska_distance.nf'
include { merge_results         } from './modules/merge_results.nf'


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
            CREATE CHANNELS FROM SAMPLESHEET
                - From the controls_ch, the samples are taxonomically classified with Sylph
                - Taxonomically classified sample reads and produces a summary of the reads
        ······································································································
        */

            Channel
                .fromPath(params.samplesheet)
                .ifEmpty { error "Sample sheet file '${params.samplesheet}' not found or empty" }
                .splitCsv(header: true, sep: ',')
                .map { row ->
                    def requiredColumns = ['sampleID', 'forward_path', 'reverse_path', 'kingdom', 'index']
                    def missingColumns = requiredColumns.findAll { !row.containsKey(it) }
                    if (missingColumns) {
                        error "Missing required column(s) in samplesheet: ${missingColumns.join(', ')}"
                    }
                        
            // Check for empty paths
                if (!row.forward_path.trim() || !row.reverse_path.trim()) {
                    error "Empty file path found for sample ${row.sampleID}. Both forward and reverse paths must be provided."
                    }
                        
                // Use the file function with error checking for the existence of the files
                def forwardFile = file(row.forward_path.trim(), checkIfExists: true)
                def reverseFile = file(row.reverse_path.trim(), checkIfExists: true)
                        
                tuple(row.sampleID.trim(), 
                    forwardFile, 
                    reverseFile, 
                    row.kingdom.trim(),
                    row.index.trim()
                    )
                }
                .branch {
                    bacteria: it[3] == 'bacteria'
                    viral: it[3] == 'viral'
                    fungal: it[3] == 'fungal'
                }
                .set { branched_samples_by_kingdom }

            // Remove the 'type' from the tuples and ensure only 3 elements
            bacteria_ch = branched_samples_by_kingdom.bacteria.map { it -> 
                tuple(it[0], it[1], it[2]) // keep only the sampleID, forward and reverse reads
            }

                bacteria_merged = bacteria_ch
                                .collect()
                                .map { samples ->
                                    def kingdom = 'bacteria'
                                    def sample_ids = samples.collect { it[0] }
                                    def forward_reads = samples.collect { it[1] }
                                    def reverse_reads = samples.collect { it[2] }
                                    tuple(sample_ids, forward_reads, reverse_reads)
                                }

            viral_ch = branched_samples_by_kingdom.viral.map { it -> 
                tuple(it[0], it[1], it[2]) // keep only the sampleID, forward and reverse reads
            }

                viral_merged = viral_ch
                                .collect()
                                .map { samples ->
                                    def kingdom = 'viral'
                                    def sample_ids = samples.collect { it[0] }
                                    def forward_reads = samples.collect { it[1] }
                                    def reverse_reads = samples.collect { it[2] }
                                    tuple(sample_ids, forward_reads, reverse_reads)
                                }

            fungal_ch = branched_samples_by_kingdom.fungal.map { it -> 
                tuple(it[0], it[1], it[2]) // keep only the sampleID, forward and reverse reads
            }

                fungal_merged = fungal_ch
                                .collect()
                                .map { samples ->
                                    def kingdom = 'fungal'
                                    def sample_ids = samples.collect { it[0] }
                                    def forward_reads = samples.collect { it[1] }
                                    def reverse_reads = samples.collect { it[2] }
                                    tuple(sample_ids, forward_reads, reverse_reads)
                                }

            bacteria_ch.view { sampleID, forward, reverse ->
                    "${color_cyan}Bacteria ctrl: ${color_green}$sampleID${color_reset} | ${color_cyan}Forward: ${color_green}$forward${color_reset} | ${color_cyan}Reverse: ${color_green}$reverse${color_reset}"
                }

            viral_ch.view { sampleID, forward, reverse ->
                    "${color_red}Viral ctrl: ${color_green}$sampleID${color_reset} | ${color_red}Forward: ${color_green}$forward${color_reset} | ${color_red}Reverse: ${color_green}$reverse${color_reset}"
                }

            fungal_ch.view { sampleID, forward, reverse ->
                    "${color_red}Fungal ctrl: ${color_green}$sampleID${color_reset} | ${color_red}Forward: ${color_green}$forward${color_reset} | ${color_red}Reverse: ${color_green}$reverse${color_reset}"
                }

        /*
        ······································································································
            NEGATIVE CONTROL WORKFLOW (NEGATIVE_CONTROL_WF)
                - From the controls_ch, the samples are taxonomically classified with Sylph
                - Taxonomically classified sample reads and produces a summary of the reads
        ······································································································
        */

        /*
        Run sylph on the reads and get read taxonomy
        */

            sylph_read_taxonomy( bacteria_merged )
            sylph_read_taxonomy( viral_merged )
            sylph_read_taxonomy( fungal_merged )

            sylph_tax_res_merge = sylph_read_taxonomy.out.collect()

        /*
        Run SKA on the reads and get kmer profiles
        */

            ska_build_reads( params.samplesheet )
            ska_distance( ska_merge.out.ska_build )
            
        /*
        Merge results and produce internative HTML to show potential contamination
        */

            merge_results( ska_distance.out.ska_merged_res, 
                            sylph_tax_res_merge )
        */




}

/*
author : Poppy J Hesketh Best
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


