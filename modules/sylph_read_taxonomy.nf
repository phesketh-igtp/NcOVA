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
        val(runID)
        file(samplesheet)

    output:
    tuple file("sylph_bacterial-merged_abundance_file.tsv"), 
        file("sylph_bacterial-merged_abundance_file.tsv"),   emit: sylph_merged_res

    script:

        """
        mkdir -p sylph_results/

        while read sampleID R1 R2; do

            echo "Processing sample: \${sampleID}"
            ln -s fastq/\${R1} \${sampleID}_R1.fastq.gz
            ln -s fastq/\${R2} \${sampleID}_R2.fastq.gz
        
        done < ${samplesheet}; do
        
        # Bacterial taxonomy
        sylph profile ${params.bacteria_sylph_db} \\
            -1 *R1*.fastq.gz -2 *R2*.fastq.gz \\
            > bacteria.tsv

        # Viral taxonomy
        sylph profile ${params.viral_sylph_db} \\
            -1 *{r,R}1*.fastq.gz -2 *{r,R}2*.fastq.gz \\
            > viral.tsv

        # Perform the taxonomic correction
        sylph-tax taxprof bacteria.tsv \\
            -t ${params.bacteria_sylph_db_id} \\
            -o sylph_results/bact-

        sylph-tax taxprof viral.tsv \\
            -t ${params.viral_sylph_db_id} \\
            -o sylph_results/vir-

        # Merge the results into a single file
                sylph-tax merge sylph/bact-*.sylphmpa \
                    --column sequence_abundance \
                    -o sylph_bact-merge_seq-abundance.tsv

                sylph-tax merge sylph/vir-*.sylphmpa \
                    --column sequence_abundance \
                    -o sylph_vir-merge_seq-abundance.tsv

        # Housekeeping of the results files
        sed -i 's/_R1.fastq.gz//g' sylph_*-merge_seq-abundance.tsv
        sed -i 's/_R2.fastq.gz//g' sylph_*-merge_seq-abundance.tsv
        sed -i 's/|/;/g' sylph_*-merge_seq-abundance.tsv

        Rscript -e ""
        library(tidyverse)

        bact <- read_tsv("sylph_bact-merge_seq-abundance.tsv", col_names = TRUE) |>
                        separate_wider_delim(clade_name, 
                        names = c("domain", "phylum", "class", 
                                "order", "family", "genome", 
                                "species", "genome"), 
                        delim = ';', too_few = 'align_start',
                        too_many = 'debug', remove = TRUE) |> 
                mutate(across(where(is.character), ~ na_if(gsub('^"|"$', '', .x), "")))
        write_tsv(bact, "sylph_bact-merge_seq-abundance.final.tsv")

        vir <- read_tsv("sylph_vir-merge_seq-abundance.tsv", col_names = TRUE)  |>
                        separate_wider_delim(clade_name, 
                        names = c("realm", "kingdom", "phylum", 
                                "class", "order", "family", 
                                "genus", "species", "genome"), 
                        delim = ';', too_few = 'align_start',
                        too_many = 'debug', remove = TRUE) |> 
                mutate(across(where(is.character), ~ na_if(gsub('^"|"$', '', .x), "")))
        write_tsv(vir, "sylph_viral-merge_seq-abundance.final.tsv")
        "
        """
}