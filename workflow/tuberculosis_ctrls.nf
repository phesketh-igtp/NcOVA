include { tbprofiler_ctrl }             from '../modules/mtbcs_controls/tbprofiler.nf'
include { mtbseq_ctrl }                 from '../modules/mtbcs_controls/mtbseq.nf'
include { tbprofiler_ctr_compile }      from '../modules/mtbcs_controls/tbprofiler_ctr_compile/main.nf'
include { mtbseq_ctr_compile }          from '../modules/mtbcs_controls/mtbseq_compile/main.nf'
include { compile_mtbc_ctrl_summary }   from '../modules/mtbcs_controls/compile_mtbc_ctrl_summary/main.nf'

workflow tuberculosis_wf {

    take:
    samplesheet

    main:

        Channel
                .fromPath(samplesheet)
                .ifEmpty { error "Sample sheet file '${samplesheet}' not found or empty" }
                .splitCsv(header: true, sep: ',')
                .map { row ->
                    def requiredColumns = ['sampleID', 'forward_path', 'reverse_path', 'species', 'index_set', 'index']
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
                        row.species.trim(),
                        row.index_set.trim(),
                        row.index.trim() // Added type for branching
                        )
                    }
                .branch {
                    mtbc: it[3] == 'control'
                }
                .set { branched_ctrls }

        /*
        ······································································································
            NEGATIVE CONTROL WORKFLOW (NEGATIVE_CONTROL_WF)
                - From the controls_ch, the samples are taxonomically classified with Sylph
                - Taxonomically classified sample reads and produces a summary of the reads
        ······································································································
        */

        // TB profiler
            tbprofiler_ctrl( branched_ctrls)
            tbprofiler_handover_collected = tbprofiler_ctrl.out.tbprofiler_results.collect()
            tbprofiler_ctr_compile( tbprofiler_handover_collected )

        // MTBSeq
            mtbseq_ctrl( branched_ctrls)

                mtbseq_class_collected = mtbseq_ctrl.out.mtbseq_class.collect()
                mtbseq_stats_collected = mtbseq_ctrl.out.mtbseq_stats.collect()

        // Compile results of MTBC control check
            compile_mtbc_ctrl_summary( mtbseq_class_collected,
                                        mtbseq_stats_collected,
                                        tbprofiler_ctr_compile.out.tbprofile_compiled
                                    )

}

/*
author: Poppy J Hesketh Best
date: 2025-04-04
version: 1.0.0-beta
description: 
    This is the tuberculosis controls workflow.
    This put the MTBC controls through the TBProfiler and MTBSeq pipelines.
changelog
        - 2024-11-01: Initial version
*/


