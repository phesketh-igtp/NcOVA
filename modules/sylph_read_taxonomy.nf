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
    
    tag "${sampleID}"

    conda params.conda_taxonomy

    publishDir "${params.outDir}/${params.runID}/sylph/", mode: 'copy'

    input:
        tuple val(sampleID), 
            path(forward), 
            path(reverse)

    output:
        file("merged_abundance_file.tsv"),   emit: sylph_merged_res

    script:

        """
        mkdir -p sylph_results/
        
        # Bacterial taxonomy
        sylph profile ${params.bacteria_sylph_db} \\
            -1 *{r,R}1*.fastq.gz -2 *{r,R}2*.fastq.gz \\
            > bacteria.tsv

        # Viral taxonomy
        sylph profile ${params.viral_sylph_db} \\
            -1 *{r,R}1*.fastq.gz -2 *{r,R}2*.fastq.gz \\
            > viral.tsv

        # Eukaryotic taxonomy
        sylph profile ${params.eukaryotic_sylph_db} \\
            -1 *{r,R}1*.fastq.gz -2 *{r,R}2*.fastq.gz \\
            > fungal.tsv

        # Perform the taxonomic correction
        sylph-tax taxprof bacteria.tsv \\
            -t ${params.bacteria_sylph_db_id} \\
            -o sylph_results/bact-

        sylph-tax taxprof viral.tsv \\
            -t ${params.viral_sylph_db_id} \\
            -o sylph_results/vir-

        sylph-tax taxprof fungal.tsv \\
            -t ${params.eukaryotic_sylph_db_id} \\
            -o sylph_results/fun-

        # Merge the results into a single file
        sylph-tax merge sylph_results/bact-*.sylphmpa \\
            --column relative_abundance \\
            -o bact-merge_abundance_file.tsv

        sylph-tax merge sylph_results/vir-*.sylphmpa \\
            --column relative_abundance \\
            -o vir-merge_abundance_file.tsv

        sylph-tax merge sylph_results/fun-*.sylphmpa \\
            --column relative_abundance \\
            -o fun-merge_abundance_file.tsv

        cat bact-merge_abundance_file.tsv \
            fun-merge_abundance_file.tsv \
            vir-merge_abundance_file.tsv \
            > merged_abundance_file.tsv

    """
}