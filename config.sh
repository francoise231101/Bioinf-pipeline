#!/bin/bash

# Configuration settings

# Log file path
LOGFILE="/home/brook/share/logfile-$(date +"%Y-%m-%d").txt"

# Directory paths
SEARCH_DIR=~/apps
SHARE_DIR=~/share
WORK_DIR=~/work

# Result directories
RESULT_QC=""
ANALYSISDIR=""
RAWDATADIR=""
BASEDIR=""
REFDIR=""
REF_FILE=""
QC_FLAG=""
TRIM_FLAG=""
INDEX_FLAG=""
MAPPING_FLAG=""
SORT_PICARD_FLAG=""
SORT_SAMTOOLS_FLAG=""
MARK_DUPLICATES_FLAG=""
TRIMMOMATIC_JAR=""
PICARD_JAR=""

# Function to display error messages and exit
exit_with_error() {
    echo "Error: $1" >&2
    exit 1
}

# Function to log messages to the log file
log_message() {
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $1" >> "$LOGFILE"
}

# Function to check if a directory exists and display an error message if not
check_directory_exists() {
    if [ ! -d "$1" ]; then
        exit_with_error "Directory not found: $1"
    fi
}

# Function to create a directory if it doesn't exist
create_directory() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
    fi
}

# Function to create a flag file
create_flag() {
    touch "$1"
}

# Function to check if a flag file exists
check_flag() {
    [ -f "$1" ]
}

# Function to check if a file has content and display an error message if not
check_file_has_content() {
    if [ ! -s "$1" ]; then
        exit_with_error "File is empty: $1"
    fi
}

# Function to set up directory paths and flags
setup_paths() {
    ANALYSISDIR="$1"
    RAWDATADIR="$2"
    BASEDIR=$(basename "$RAWDATADIR")

    # Flag files for each step
    QC_FLAG="$RESULT_QC/qc_completed.flag"
    TRIM_FLAG="$RESULT_QC/trim_completed.flag"
    INDEX_FLAG="$SHARE_DIR/index_completed.flag"
    MAPPING_FLAG="$RESULT_QC/mapping_completed.flag"
    SORT_PICARD_FLAG="$RESULT_QC/sort_picard_completed.flag"
    SORT_SAMTOOLS_FLAG="$RESULT_QC/sort_samtools_completed.flag"
    MARK_DUPLICATES_FLAG="$RESULT_QC/duplication_removed_completed.flag"

    RESULT_QC="$ANALYSISDIR/$BASEDIR"
    create_directory "$RESULT_QC"
    create_directory "$RESULT_QC/RESULT_TRIM"
    create_directory "$RESULT_QC/RESULT_MAPPING"

    REFDIR="$SHARE_DIR/ref"
    create_directory "$REFDIR"
    REF_FILE="$REFDIR/hg38.fa"
    PREFIX_REF="$REFDIR/hg38_BWA_idx"

    # Check if reference file or directory exists, skip if necessary
    if [ ! -e "$REF_FILE" ] && [ ! -d "$REFDIR" ]; then
        log_message "Downloading reference file"
        wget -P "$REFDIR" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz || exit_with_error "Failed to download reference file."

        # Check for incomplete download or corrupted file
        if [ -s "$REFDIR/hg38.fa.gz" ]; then
            gunzip -c "$REFDIR/hg38.fa.gz" > "$REFDIR/hg38.fa" || exit_with_error "Failed to unzip reference file."
        else
            exit_with_error "Downloaded file is empty or corrupted. Please try downloading again."
        fi
    else
        log_message "Reference file or directory already exists"
    fi

    # Check if BWA index flag is not set
    if ! check_flag "$INDEX_FLAG"; then
        log_message "Starting Indexing of Human genome reference"

        # Check if any file with the specified prefix exists
        file_count=$(find "$REFDIR" -maxdepth 1 -type f -name "$(basename "$PREFIX_REF").*" | wc -l)

        if [ "$file_count" -gt 0 ]; then
            log_message "Skipping indexing as it has been completed before"
        else
            bwa index -p "$PREFIX_REF" "$REFDIR/hg38.fa" || exit_with_error "Failed to create BWA index."
            create_flag "$INDEX_FLAG"
        fi
    else
        log_message "Skipping indexing as it has been completed before"
    fi

    # Trimmomatic JAR
    TRIMMOMATIC_JAR=$(find "$SEARCH_DIR" -name "trimmomatic*.jar" -type f -print -quit)
    [ -n "$TRIMMOMATIC_JAR" ] || exit_with_error "Trimmomatic JAR not found in $SEARCH_DIR or its subdirectories."

    # Picard JAR
    PICARD_JAR=$(find "$SEARCH_DIR" -name "picard*.jar" -type f -print -quit)
    [ -n "$PICARD_JAR" ] || exit_with_error "Picard JAR not found in $SEARCH_DIR or its subdirectories."
}
