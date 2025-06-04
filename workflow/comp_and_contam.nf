include { sylph_read_taxonomy   } from '../modules/sylph_read_taxonomy.nf'
include { seqkit_stats          } from '../modules/seqkit_stats.nf'
include { ska_reads_build       } from '../modules/read_ska.nf'
include { ska_distance          } from '../modules/ska_distance.nf'
include { merge_results         } from '../modules/merge_results.nf'

workflow comp_and_contam_wf {

    take:
    samplesheet

    main:

    //def color_purple = '\u001B[35m'
    def color_green  = '\u001B[32m'
    def color_red    = '\u001B[31m'
    def color_reset  = '\u001B[0m'
    def color_cyan   = '\u001B[36m'

            Channel
                .fromPath(samplesheet)
                .ifEmpty { error "Sample sheet file '${samplesheet}' not found or empty" }
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
        Run sylph on the reads and get read taxonomy
        */
            sylph_read_taxonomy( samplesheet )
            sylph_tax_res_merge = sylph_read_taxonomy.out.collect()

        /*
        Run SKA on the reads and get kmer profiles
        */
            ska_reads_build( samplesheet )
            ska_distance( ska_reads_build.out.ska_build )

        /*
        Merge results and produce internative HTML to show potential contamination
        */
            merge_results( ska_distance.out.ska_merged_res, 
                            sylph_tax_res_merge )

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


