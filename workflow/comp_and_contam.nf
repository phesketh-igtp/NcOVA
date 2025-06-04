include { sylph_read_taxonomy   } from '../modules/sylph_read_taxonomy.nf'
include { seqkit_stats          } from '../modules/seqkit_stats.nf'
include { ska_reads_build       } from '../modules/read_ska.nf'
include { ska_distance          } from '../modules/ska_distance.nf'
include { merge_results         } from '../modules/merge_results.nf'

workflow comp_and_contam_wf {

    take:
    samplesheet

    main:

        /*
        Run sylph on the reads and get read taxonomy
        */
            sylph_read_taxonomy( samplesheet )

        /*
        Run SKA on the reads and get kmer profiles
        */
            ska_reads_build( samplesheet )
            ska_distance( ska_reads_build.out.ska_build )

        /*
        Merge results and produce internative HTML to show potential contamination
        */
            merge_results( ska_distance.out.ska_distance, 
                            sylph_read_taxonomy.out.sylph_merged_res )

}

/*
author: Poppy J Hesketh Best
date: 2025-06-04
version: 1.0.0-beta
description: 
changelog:
        v1.0.0-2025-06-04: Initial version
*/


