process sylph_classification {

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
        path(sylph_profiles)

    output:
    path("merged_sylph_sequence_abundance_file.tsv"),   emit: sylph_merged_res

    script:

        """
        # Get taxonomy for the profiles
            sylph-tax download --download-to my_existing_folder/

        # Classify using the Sylph taxonomy database   
            sylph-tax taxprof sylph_profiles \\
                -t ${params.sylph_bbdd_id} \\
                -o sylph/tax_

        # remove any empty files
            find sylph/ -type f -name 'tax_*.sylphmpa' -empty -delete

        # Merge the results
            sylph-tax merge \\
                sylph/tax_*.sylphmpa \\
                --column sequence_abundance \\
                -o merged_sylph_sequence_abundance_file.tsv
        """
}