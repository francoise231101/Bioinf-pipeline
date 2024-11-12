#!/bin/bash

# Check if ANALYSISDIR and case_name are passed as arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: sh $0 /path/to/analysisdir case_name"
    exit 1
fi

# Set ANALYSISDIR and case_name from the arguments
ANALYSISDIR="$1"
case_name="$2"

# Verify that the ANALYSISDIR exists
if [ ! -d "$ANALYSISDIR" ]; then
    echo "Error: Directory '$ANALYSISDIR' does not exist."
    exit 1
fi

# Generate the output directory name with case name and timestamp
timestamp=$(date +"%Y%m%d-%H:%M:%S")
output_dir="${case_name}-${timestamp}"

# Create the output directory
mkdir -p "$output_dir"

# Set log file in the output directory
LOGFILE="$output_dir/preprocess_wgs.$timestamp.log"
exec 1>> "$LOGFILE" 2>&1


# Sample list
samples=$(find "$ANALYSISDIR" -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)

# Organize the folder & rename files
echo "### Organizing directory & setting analysis parameters ###"
echo "Start: $(date)"
echo "End: $(date)"