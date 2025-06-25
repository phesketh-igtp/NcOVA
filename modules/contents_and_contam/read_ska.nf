process ska_reads_build {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-04
    @version: 1.0.2
    @description: 
        Run Sylph on the reads and get read taxonomy and statistics of the reads. Outputs of this modules 
        are intended to be combines into a single file for each sample, and then concatenated into a
        single file for all samples for a particular run.
    @chagelog: 
        v1.0.2-2025-04-04: created new paths for the results `Classification/ Statistics/`
        v2.0.0-2025-05-27: updated to use Sylph instead of Kraken2, and added support for empty controls
*/
    
    conda params.general_env

    publishDir "${params.outDir}/${params.runID}/ska/", mode: 'copy'

    input:
        path(samplesheet)

    output:
        path("seqs.skf"), emit: ska_build

    script:

        """
        sed 's/,/\t/g' ${samplesheet} \\
            | sed '1d' \\
            | cut -f1,2,3 > input_sequence.txt

        line_count=\$(wc -l < input_sequence.txt)

        if [ "\$line_count" -le 19 ]; then threads=1
            elif [ "\$line_count" -le 39 ]; then threads=2
            elif [ "\$line_count" -le 79 ]; then threads=4
            else threads=8
        fi; echo "Using \$threads threads"

        wc -l input_sequence.txt | \\
            awk '{print \$1-1}' > num_samples.txt

        ska build -f input_sequence.txt \\
            --min-count 5 \\
            --min-qual 20 \\
            --qual-filter strict \\
            --threads \$threads \\
            -o seqs
        """
}