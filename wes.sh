#!/bin/bash

### AsiaMed Common Pipeline ###                       
## commun_wes.sh version 0.1
## Description: for exome data
## Usage: sh ~apps/commun_wes.sh ANALYSISDIR=/path/to/analysisdir RAWDATADIR=/path/to/rawdatadir
## Output: no standard output
## Requirements: Require pipeline scripts    

## Author: Francoise (Kuan-Hua) TU ART.
## Creation Date:
## Last revision date:
## Known bugs:

# Include the configuration file
. /home/francoise/apps/config.sh

# Redirect stdout and stderr to the log file
exec &> "$LOGFILE"

# Log the start of the pipeline
log_message "Starting the pipeline"

# Check command-line arguments
if [ "$#" -ne 2 ]; then
    exit_with_error "Usage: $0 '/path/to/analysisdir/' '/path/to/rawdatadir/'"
fi

ANALYSISDIR="$1"
RAWDATADIR="$2"

check_directory_exists "$RAWDATADIR"

# Set up paths and flags
setup_paths "$ANALYSISDIR" "$RAWDATADIR"

#########################################################################################
## Step 1. QC (Quality Control)                                                      ##
#########################################################################################

# Run fastqc only if the QC step has not been completed before
if ! check_flag "$QC_FLAG"; then
    # Log QC step
    log_message "Starting QC step"
    log_message "run fastqc *.fastq.gz -o "$RESULT_QC" "

    # Run fastqc on each .fastq.gz file in RAWDATADIR using parallel
    find "$RAWDATADIR" -name "*.fastq.gz" | parallel 'fastqc -o '"$RESULT_QC"' {}'

    # Log QC completion
    log_message "QC step completed"

    # Create a flag file to indicate that QC is completed
    create_flag "$QC_FLAG"
else
    log_message "Skipping QC step as it has been completed before"
fi

#########################################################################################
## Step 2. Trimming by using Trimmomatic                                               ##
#########################################################################################

# Run trimmomatic only if the trimming step has not been completed before
if ! check_flag "$TRIM_FLAG"; then
    # Log Trimming step
    log_message "Starting Trimming step"

    # Run trimmomatic on each *fastq.gz file in RAWDATADIR without parallel
    for r1_file in "$RAWDATADIR"/*_R1.fastq.gz; do
        n=$(basename "$r1_file" "_R1.fastq.gz")
        r2_file="$RAWDATADIR/${n}_R2.fastq.gz"
        log_file="$RESULT_QC/RESULT_TRIM/${n}_trimmed.log"

        log_message "java -jar "$TRIMMOMATIC_JAR" PE -threads 15 -phred64 \
                        "$r1_file" "$r2_file" \
                        "$RESULT_QC/RESULT_TRIM/${n}_R1_trimmed.fastq.gz" \
                        "$RESULT_QC/RESULT_TRIM/${n}_R1_unpaired.fastq.gz" \
                        "$RESULT_QC/RESULT_TRIM/${n}_R2_trimmed.fastq.gz" \
                        "$RESULT_QC/RESULT_TRIM/${n}_R2_unpaired.fastq.gz" \
                        ILLUMINACLIP:"$SEARCH_DIR"/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
                        SLIDINGWINDOW:4:20 MINLEN:15 "

        if [ -f "$r2_file" ]; then
            log_message "Processing $n"
            java -jar "$TRIMMOMATIC_JAR" PE -threads 15 -phred64 \
                "$r1_file" "$r2_file" \
                "$RESULT_QC/RESULT_TRIM/${n}_R1_trimmed.fastq.gz" \
                "$RESULT_QC/RESULT_TRIM/${n}_R1_unpaired.fastq.gz" \
                "$RESULT_QC/RESULT_TRIM/${n}_R2_trimmed.fastq.gz" \
                "$RESULT_QC/RESULT_TRIM/${n}_R2_unpaired.fastq.gz" \
                ILLUMINACLIP:"$SEARCH_DIR"/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
                SLIDINGWINDOW:4:20 MINLEN:15 &> "$log_file"
        else
            log_message "Error: Missing corresponding R2 file for $r1_file"
        fi
    done

    # Log Trimming completion
    log_message "Trimming step completed"

    # Create a flag file to indicate that trimming is completed
    create_flag "$TRIM_FLAG"
else
    log_message "Skipping Trimming step as it has been completed before"
fi

#########################################################################################
## Step 3: Indexing Human genome reference                                             ##
#########################################################################################
# PREFIX_REF="$REFDIR/hg38_BWA_idx"

# if ! check_flag "$INDEX_FLAG"; then
#     log_message "Starting Indexing of Human genome reference"

#     if [ ! -f "$PREFIX_REF".sa ]; then
#         bwa index -p "$PREFIX_REF" "$REF_FILE" || exit_with_error "Failed to create BWA index."
#     else
#         log_message "Skipping indexing as it has been completed before"
#         create_flag "$INDEX_FLAG"
#     fi
# else
#     log_message "Skipping indexing as it has been completed before"
# fi

#########################################################################################
## Step 4. Mapping to human genome reference GRCh38 using BWA                          ##
#########################################################################################

if ! check_flag "$MAPPING_FLAG"; then
    # Mapping the reads
    log_message "Starting mapping step"

    for r1_trimmed_file in "$RESULT_QC"/RESULT_TRIM/*_R1_trimmed.fastq.gz; do
        n=$(basename "$r1_trimmed_file" "_R1_trimmed.fastq.gz")
        r2_trimmed_file="$RESULT_QC/RESULT_TRIM/${n}_R2_trimmed.fastq.gz"
        mapping_sam_file="$RESULT_QC/RESULT_MAPPING/${n}_mapping.sam"
        log_file="$RESULT_QC/RESULT_MAPPING/${n}_mapping.log"

        log_message "bwa mem -t 20 "$PREFIX_REF" "$RESULT_QC/RESULT_TRIM/${n}_R1_trimmed.fastq.gz" "$RESULT_QC/RESULT_TRIM/${n}_R2_trimmed.fastq.gz" > "$mapping_sam_file" "

        if [ -f "$r2_trimmed_file" ]; then
            log_message "Processing $n"
            bwa mem -t 20 "$PREFIX_REF" \
                "$RESULT_QC/RESULT_TRIM/${n}_R1_trimmed.fastq.gz" "$RESULT_QC/RESULT_TRIM/${n}_R2_trimmed.fastq.gz" > "$mapping_sam_file" &> "$log_file"
        else
            log_message "Error: Missing corresponding R2 trimmed file for $r1_trimmed_file"
        fi
    done

    # Log Mapping completion
    log_message "Mapping step completed"
    create_flag "$MAPPING_FLAG"

else
    log_message "Skipping mapping as it has been completed before"
fi

#########################################################################################
## Step 5 (Option 1): After mapping, sorting by picard and merging                     ##
#########################################################################################
if ! check_flag "$SORT_PICARD_FLAG"; then
    log_message "Sorting by Picard"

    # Loop through each SAM file
    for sam_file in "$RESULT_QC/RESULT_MAPPING/"*_mapping.sam; do
        n=$(basename "$sam_file" "_mapping.sam")
        log_file="$RESULT_QC/RESULT_MAPPING/${n}_sorted.log"

        log_message "java -jar SortSam INPUT="$sam_file" \
            OUTPUT="$RESULT_QC/RESULT_MAPPING/${n}_mapping_picard_sorted.bam" \
            SORT_ORDER=coordinate \
            VALIDATION_STRINGENCY=STRICT || exit_with_error "Picard SortSam failed." "

        # Sort the SAM file using Picard
        java -jar "$PICARD_JAR" SortSam \
            INPUT="$sam_file" \
            OUTPUT="$RESULT_QC/RESULT_MAPPING/${n}_mapping_picard_sorted.bam" \
            SORT_ORDER=coordinate \
            VALIDATION_STRINGENCY=STRICT || exit_with_error "Picard SortSam failed." &> "$log_file" 
    done

    log_message "Sorting completed"
    create_flag "$SORT_PICARD_FLAG"

else
    log_message "Skipping sorting as it has been completed before"
fi

#########################################################################################
## Step 5 (Option 2): After mapping, sorting by samtools                               ##
#########################################################################################
if ! check_flag "$SORT_SAMTOOLS_FLAG"; then
    log_message "Sorting by samtools sort"

    # Loop through each SAM file
    for sam_file in "$RESULT_QC/RESULT_MAPPING/"*_mapping.sam; do
        n=$(basename "$sam_file" "_mapping.sam")
        log_file="$RESULT_QC/RESULT_MAPPING/${n}_sorted.log"

        log_message "samtools view -S -b "$sam_file" -o "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools.bam""
        log_message "samtools sort -o "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools_sorted.bam" "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools.bam""
        log_message "samtools index "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools_sorted.bam""

        # Sort the SAM file using samtools sort
        samtools view -S -b "$sam_file" -o "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools.bam"
        samtools sort -o "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools_sorted.bam" "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools.bam"
        samtools index "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools_sorted.bam"
        samtools view "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools_sorted.bam" -o "$RESULT_QC/RESULT_MAPPING/${n}_mapping_samtools_sorted.sam" || exit_with_error "samtools sort failed." &> "$log_file" 
    done

    log_message "Sorting completed"
    create_flag "$SORT_SAMTOOLS_FLAG"

else
    log_message "Skipping sorting as it has been completed before"
fi

#########################################################################################
## Step 6: Picard MarkDuplicates to remove the duplications                           ##
#########################################################################################
if ! check_flag "$MARK_DUPLICATES_FLAG"; then
    log_message "Marking duplicates with Picard"

    # Loop through each sorted BAM file
    for sorted_bam in "$RESULT_QC/RESULT_MAPPING/"*_sorted.bam; do
        n=$(basename "$sorted_bam" "_sorted.bam")
        log_file="$RESULT_QC/RESULT_MAPPING/${n}_sorted_dedup.log"

        log_message "java -jar "$PICARD_JAR" MarkDuplicates \
                        CREATE_INDEX=true \
                        INPUT="$sorted_bam" \
                        OUTPUT="$RESULT_QC/RESULT_MAPPING/${n}_sorted_dedup.bam" \
                        METRICS_FILE="$RESULT_QC/RESULT_MAPPING/${n}_sorted_dedup_metrics.txt" \
                        VALIDATION_STRINGENCY=STRICT || exit_with_error "MarkDuplicates failed.""

        # Mark duplicates using Picard
        java -jar "$PICARD_JAR" MarkDuplicates \
            CREATE_INDEX=true \
            INPUT="$sorted_bam" \
            OUTPUT="$RESULT_QC/RESULT_MAPPING/${n}_sorted_dedup.bam" \
            METRICS_FILE="$RESULT_QC/RESULT_MAPPING/${n}_sorted_dedup_metrics.txt" \
            VALIDATION_STRINGENCY=STRICT || exit_with_error "MarkDuplicates failed." &> "$log_file" 
    done

    log_message "Marking duplicates completed"
    create_flag "$MARK_DUPLICATES_FLAG"

else
    log_message "Skipping MarkDuplicates as it has been completed before"
fi

#########################################################################################
## Step 7: Variant calling by MuTect2                                                 ##
#########################################################################################

if ! check_flag "$MUTECT2_FLAG"; then
    log_message "Starting MuTect2 variant calling"

    # Loop through each sorted and deduplicated BAM file
    for dedup_bam in "$RESULT_QC/RESULT_MAPPING/"*_sorted_dedup.bam; do
        n=$(basename "$dedup_bam" "_sorted_dedup.bam")
        output_vcf="$RESULT_QC/RESULT_MUTECT2/${n}_variant_calls.vcf.gz"
        log_file="$RESULT_QC/RESULT_MUTECT2/${n}_mutect2.log"

        log_message "gatk Mutect2 \
                        -R "$REFERENCE_FA" \
                        -I "$dedup_bam" \
                        -O "$output_vcf" &> "$log_file" "

        # Run MuTect2
        gatk Mutect2 \
            -R "$REFERENCE_FA" \
            -I "$dedup_bam" \
            -O "$output_vcf" &> "$log_file"

        if [ $? -eq 0 ]; then
            log_message "MuTect2 variant calling for $n completed successfully."
            # Add additional checks or actions here if needed
        else
            exit_with_error "Error in MuTect2 variant calling for $n. Check $log_file for details."
        fi
    done

    log_message "MuTect2 variant calling completed"
    create_flag "$MUTECT2_FLAG"

else
    log_message "Skipping MuTect2 variant calling as it has been completed before"
fi

# Log the end of the pipeline
log_message "Pipeline completed"
