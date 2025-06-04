process tbprofiler_ctrl {

        tag "${sampleID}"
        
        conda params.tbprofiler_env

        input:
                tuple val(sampleID), 
                        path(forward),
                        path(reverse)
                        
        output:
                tuple val(sampleID), 
						path("results/tbdb-${sampleID}.results.txt"), emit: tbprofiler_results
                        path("results/tbdb-${sampleID}.results.json"), optional: true, emit: tbprofiler_results
				
				path("log/${sampleID}_tb_profiler_status.txt"), emit: tbprofiler_handover

        script:
        """
        # Run TB-Profiler using TBDB database
            set +e #tells the shell not to exit immediately if a command fails
            tb-profiler profile \\
                    -1 ${forward} \\
                    -2 ${reverse} \\
                    -p tbdb-${sampleID} \\
                    --txt --dir . \\
                    --db ${params.outDir}/db/tbprofiler/tbdb \\
                    --threads ${task.cpus} \\
                    --no_delly \\
                    > tb-profiler.out 2> tb-profiler.err
            set -e #restores default command fails checks

        mkdir -p log/

        # Check if the results file was created
			if [[ -f results/tbdb-${sampleID}.results.txt ]]; then
				echo "${sampleID},SUCCESS" > log/${sampleID}_tb_profiler_status.txt
            else
                    echo "${sampleID},FAILED" > log/${sampleID}_tb_profiler_status.txt
                    # Create empty files to satisfy output requirements
                    echo "" > results/tbdb-${sampleID}.results.txt
                    echo "" > results/tbdb-${sampleID}.results.json
            fi

        # Always exit with status 0 to prevent pipeline failure
                exit 0
        """

}