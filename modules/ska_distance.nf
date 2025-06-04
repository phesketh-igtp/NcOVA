process ska_distance {

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
        path(ska_build)

    output:
        file("distances.txt"), emit: ska_distance
        file("clusters.txt"), emit: ska_clusters

    script:

        """
        ska distance -i ${ska_build} \\
            -o distances.txt

        python scripts/cluster_dists.py distances.txt \\
            --snps 21 \\
            --mismatches 0.05
        """
}