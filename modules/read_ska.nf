process read_ska {

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
    
    conda params.conda_taxonomy

    publishDir "${params.outDir}/${params.runID}/ska/", mode: 'copy'

    input:
        val(runID)
        file(samplesheet)

    output:
    tuple file("distances.txt"), 
        file("lo_out"),   emit: ska_merged_res

    script:

        """
        sed 's/,/\t/g' ${samplesheet} | \\
            cut -f2,3,4 > input_sequence.txt

        ska build -f input_sequence.txt \\
            --min-count 5 \\
            --min-qual 20 \\
            --qual-filter strict \\
            --threads 10 \\
            -o seqs

        ska lo seqs.skf lo_out

        ska distance -o distances.txt seqs.skf

        python scripts/cluster_dists.py \\
            distances.txt \\
            --snps 30 \\
            --mismatches 0.05
        """
}