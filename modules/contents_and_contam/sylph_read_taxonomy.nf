process sylph_read_taxonomy {

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

    publishDir "${params.outDir}/${params.runID}/sylph/", mode: 'copy'

    input:
        path(samplesheet)

    output:
    path("sylph_profile.tsv"),   emit: sylph_profiles

    script:

        """
    # Create directory
        mkdir -p sylph
    
    # Parse the reads in the samplesheet
        while IFS=',' read -r sampleID forward reverse species index; do

            ln -s "\${forward}" "\${sampleID}_R1.fastq.gz"
            ln -s "\${reverse}" "\${sampleID}_R2.fastq.gz"
                
            sylph sketch \\
                -1 "\${sampleID}_R1.fastq.gz" \\
                -2 "\${sampleID}_R2.fastq.gz" \\
                -d sylph/ \\
                -t ${task.cpus}

        done < <(sed '1d' ${params.samplesheet})

    # Profile the sketches
        sylph profile \\
                ${params.sylph_bbdd} \\
                sylph/*syl* \\
                --estimate-unknown \\
                --read-seq-id 0.99 \\
                -t ${task.cpus} \\
                -o sylph_profile.tsv
        """
}