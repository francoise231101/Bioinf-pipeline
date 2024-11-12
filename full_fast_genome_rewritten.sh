
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
# Replace 'rename_fastq.py' with actual command if needed
python /path/to/rename_fastq.py -d "$ANALYSISDIR" -f wes
python /path/to/organize_data_folder.py -d "$ANALYSISDIR" -p /path/to/pipelinebase -t /path/to/targetlist
echo "End: $(date)"

sleep 5

# Quality Control (FastQC)
echo "### First pass FastQC ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "FP fastqc for sample: $currentSample"
    # R1
    fastqc -o "$ANALYSISDIR/$currentSample/QC" "$ANALYSISDIR/$currentSample/$currentSample.R1.fastq.gz" \
        > "$ANALYSISDIR/$currentSample/logs/process_fastqc_R1.$logbasename.log" 2>&1
    # R2
    fastqc -o "$ANALYSISDIR/$currentSample/QC" "$ANALYSISDIR/$currentSample/$currentSample.R2.fastq.gz" \
        > "$ANALYSISDIR/$currentSample/logs/process_fastqc_R2.$logbasename.log" 2>&1
done
echo "End: $(date)"

# Trim FastQ Files
echo "### Trim FastQ ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "Trimming sample: $currentSample"
    # Ensure that paths to `trim_galore` or any other trimming tools are correctly set up
    trim_galore --paired "$ANALYSISDIR/$currentSample/$currentSample.R1.fastq.gz" \
                         "$ANALYSISDIR/$currentSample/$currentSample.R2.fastq.gz" \
                         --output_dir "$ANALYSISDIR/$currentSample" \
                         > "$ANALYSISDIR/$currentSample/logs/trim_fastq.$logbasename.log" 2>&1
done
echo "End: $(date)"



# Second Pass FastQC
echo "### Second pass FastQC ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "SP FastQC for sample: $currentSample"
    # R1
    fastqc -o "$ANALYSISDIR/$currentSample/QC/$currentSample.trimmed_fastqc" \
           "$ANALYSISDIR/$currentSample/$currentSample.R1.fastq.gz" \
           > "$ANALYSISDIR/$currentSample/logs/process_fastqc_SP_R1.$logbasename.log" 2>&1
    # R2
    fastqc -o "$ANALYSISDIR/$currentSample/QC/$currentSample.trimmed_fastqc" \
           "$ANALYSISDIR/$currentSample/$currentSample.R2.fastq.gz" \
           > "$ANALYSISDIR/$currentSample/logs/process_fastqc_SP_R2.$logbasename.log" 2>&1
done
echo "End: $(date)"

# Alignment using BWA-MEM
echo "### Computing read group and aligning sequences ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "Aligning sample: $currentSample"
    bwa mem -t 16 -R "@RG\tID:$currentSample\tSM:$currentSample\tPL:illumina" \
            /path/to/reference.fasta \
            "$ANALYSISDIR/$currentSample/$currentSample.R1.fastq.gz" \
            "$ANALYSISDIR/$currentSample/$currentSample.R2.fastq.gz" \
            | samtools sort -o "$ANALYSISDIR/$currentSample/$currentSample.sort.bam"
done
echo "End: $(date)"

# Clean BAM Files
echo "### Cleaning BAM files ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "Cleaning BAM for sample: $currentSample"
    samtools view -F 4 -b "$ANALYSISDIR/$currentSample/$currentSample.sort.bam" \
        > "$ANALYSISDIR/$currentSample/$currentSample.clean.bam"
done
echo "End: $(date)"

# Mark Duplicates
echo "### Mark duplicates in BAM files ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "Marking duplicates for sample: $currentSample"
    samtools markdup "$ANALYSISDIR/$currentSample/$currentSample.clean.bam" \
                     "$ANALYSISDIR/$currentSample/$currentSample.bam"
    # Store metrics file in QC directory
    echo "Mark duplicates metrics for $currentSample" > "$ANALYSISDIR/$currentSample/QC/$currentSample.dedup.metrics"
done
echo "End: $(date)"

# Check BAM Files
echo "### Checking BAM files ###"
echo "Start: $(date)"
for currentSample in $samples; do
    echo "Checking BAM integrity for sample: $currentSample"
    samtools quickcheck "$ANALYSISDIR/$currentSample/$currentSample.bam" || \
        echo "$currentSample failed BAM integrity check" >> "$ANALYSISDIR/$currentSample/defective_bam.list"
done
echo "End: $(date)"

# Configure Variant Calling (Sample List Preparation)
echo "### Configuring for Variant Calling ###"
echo "Start: $(date)"

# Control sample list
if [ -e "$ANALYSISDIR/control.list" ]; then
    rm "$ANALYSISDIR/control.list"
fi
# Generate control list manually if available or skip if not applicable

# Batch sample list
if [ -e "$ANALYSISDIR/batch.list" ]; then
    rm "$ANALYSISDIR/batch.list"
fi
find "$ANALYSISDIR" -maxdepth 1 -mindepth 1 -type d -exec basename {} \; >> "$ANALYSISDIR/batch.list"

# List of samples not in control dataset
if [ -e "$ANALYSISDIR/gvcf.list" ]; then
    rm "$ANALYSISDIR/gvcf.list"
fi

for currentSample in $samples; do
    if grep -q "$currentSample" "$ANALYSISDIR/control.list"; then
        echo "$currentSample is part of control dataset"
    else
        echo "$ANALYSISDIR/$currentSample/$currentSample.raw.g.vcf.gz" >> "$ANALYSISDIR/gvcf.list"
    fi
done

echo "End: $(date)"


# Configure variant calling analysis
echo "### Configuration step for Variant Calling ###"
echo "Start: $(date +"%F_%H-%M-%S")"

# Clear previous lists if they exist
varcallfiles=("$ANALYSISDIR"/*/*.varcall.list)
if [ -e "${varcallfiles[0]}" ]; then
    rm "$ANALYSISDIR"/*/*.varcall.list
fi

denovofiles=("$ANALYSISDIR"/*/*.denovo)
if [ -e "${denovofiles[0]}" ]; then
    rm "$ANALYSISDIR"/*/*.denovo
fi

# Check if family files exist
situation=0
familylist=("$ANALYSISDIR"/*.family.list)
if [ -e "${familylist[0]}" ]; then
    situation=2
fi

# If no family file, handle each sample individually
if [ "$situation" -eq "0" ]; then
    echo "No family file found. Processing each sample independently."
    for currentSample in $samples; do
        echo "$currentSample" >> "$ANALYSISDIR/$currentSample/$currentSample.varcall.list"
    done
fi

# If family files exist, handle accordingly
if [ "$situation" -eq "2" ]; then
    echo "Family files found. Processing family samples."

    # Collect all individuals from family files
    declare -a all_indiv_family_tab
    for familyfile in "$ANALYSISDIR"/*.family.list; do
        while IFS= read -r line; do
            all_indiv_family_tab+=("$line")
        done < "$familyfile"
    done

    # Process each sample individually if not in any family file
    for currentSample in $samples; do
        found=0
        for a in "${all_indiv_family_tab[@]}"; do
            if [ "$a" == "$currentSample" ]; then
                found=1
                break
            fi
        done
        if [ "$found" -eq 0 ]; then
            echo "$currentSample" >> "$ANALYSISDIR/$currentSample/$currentSample.varcall.list"
            echo "Individual sample: $currentSample"
        fi
    done

    # Process each family
    for familyfile in "$ANALYSISDIR"/*.family.list; do
        echo "Processing family file: $familyfile"

        # Initialize family members
        unset family_tab casindex ci mere pere
        while IFS= read -r line; do
            # Check relation type
            parente=$(grep "Parente" "$ANALYSISDIR/$line/$line.info" | cut -f2)
            if [ "$parente" == "Cas index" ]; then
                ci="$line"
            elif [ "$parente" == "Mère" ]; then
                mere="$line"
            elif [ "$parente" == "Père" ]; then
                pere="$line"
            fi
            family_tab+=("$line")
        done < "$familyfile"

        # If all family members (case, mother, father) are found, set for de novo analysis
        if [ -n "$ci" ] && [ -n "$mere" ] && [ -n "$pere" ]; then
            echo "$ci" >> "$ANALYSISDIR/$ci/$ci.denovo"
        fi

        # Process family for variant calling
        casindex_tab=("${family_tab[@]}")
        for casindex in "${casindex_tab[@]}"; do
            # Add case index sample to variant call list
            echo "$casindex" >> "$ANALYSISDIR/$casindex/$casindex.varcall.list"
            echo "Cas index sample: $casindex"

            # Add other family members to the case index's variant call list
            for member in "${family_tab[@]}"; do
                if [ "$member" != "$casindex" ]; then
                    echo "$member" >> "$ANALYSISDIR/$casindex/$casindex.varcall.list"
                    echo "Family sample: $member"
                fi
            done
        done
    done
fi

echo "End: $(date +"%F_%H-%M-%S")"


echo "### Final Report Compilation Step ###"
echo "Start: $(date +"%F_%H-%M-%S")"

# Set up the final report output location
REPORT_DIR="$ANALYSISDIR/final_report"
mkdir -p "$REPORT_DIR"

# Compile CNV plots and reports for each sample
for currentSample in $samples; do
    echo "Compiling final report for sample: $currentSample"
    
    # Copy CNV annotations and plots to the report directory
    cp "$ANALYSISDIR/$currentSample/$currentSample.cnv.annot.tsv" "$REPORT_DIR/${currentSample}_cnv_annotations.tsv"
    cp "$ANALYSISDIR/$currentSample/$currentSample.cnv.annot.plot.png" "$REPORT_DIR/${currentSample}_cnv_plot.png"
    cp "$ANALYSISDIR/$currentSample/$currentSample.lumpy.vcf" "$REPORT_DIR/${currentSample}_lumpy.vcf"
    cp "$ANALYSISDIR/$currentSample/$currentSample.controlfreec.tsv" "$REPORT_DIR/${currentSample}_controlfreec.tsv"

    echo "Successfully compiled CNV plots and reports for $currentSample"
done

# Aggregate all CNV annotations into a single report
echo "Generating consolidated CNV report..."
cat "$REPORT_DIR"/*_cnv_annotations.tsv > "$REPORT_DIR/consolidated_cnv_annotations.tsv"
echo "Consolidated CNV report generated: $REPORT_DIR/consolidated_cnv_annotations.tsv"

# Final cleanup and status output
echo "Final report generation complete."
echo "Pipeline processing finished for all samples."
echo "End: $(date +"%F_%H-%M-%S")"


# Variant Annotation Step
echo "### Variant Annotation ###"
echo "Start: $(date +"%F_%H-%M-%S")"
for varcallfile in $(ls "$ANALYSISDIR"/*/*.varcall.list); do
    casindex_sample=$(basename "${varcallfile%%.*}")
    reportbasename=""

    # Check and clean previous annotation files
    if [ -e "$ANALYSISDIR/$casindex_sample/$casindex_sample.*.annot.vcf" ]; then
        rm "$ANALYSISDIR/$casindex_sample/$casindex_sample.*.annot.vcf"
    fi

    # Run the variant annotation process
    INPUTFILE="$ANALYSISDIR/$casindex_sample/$casindex_sample.snpeff.vcf"
    OUTPUTFILE="$ANALYSISDIR/$casindex_sample/$casindex_sample.annot.vcf"
    LOGFILE="$ANALYSISDIR/$casindex_sample/logs/annotate_variants.$logbasename.log"
    CNVREPORT="$ANALYSISDIR/$casindex_sample/$casindex_sample.cnv.annot.tsv"

    echo "Annotating variants for sample: $casindex_sample"
    bash "$PIPELINEBASE/common/vcf/wrapper_annotate_variants.sh" \
        -i "$INPUTFILE" \
        -o "$OUTPUTFILE" \
        -c "$CNVREPORT" \
        &>> "$LOGFILE"
done
echo "End: $(date +"%F_%H-%M-%S")"

# Rare Variant Analysis Step
echo "### Rare Variants Analysis ###"
echo "Start: $(date +"%F_%H-%M-%S")"
for varcallfile in $(ls "$ANALYSISDIR"/*/*.varcall.list); do
    casindex_sample=$(basename "${varcallfile%%.*}")

    # Clean previous rare variant files if they exist
    if [ -e "$ANALYSISDIR/$casindex_sample/$casindex_sample.*.rare.vcf" ]; then
        rm "$ANALYSISDIR/$casindex_sample/$casindex_sample.*.rare.vcf"
    fi

    # Prepare inputs for rare variant analysis
    INPUTFILE="$ANALYSISDIR/$casindex_sample/$casindex_sample.annot.vcf"
    OUTPUTFILE="$ANALYSISDIR/$casindex_sample/$casindex_sample.rare.vcf"
    LOGFILE="$ANALYSISDIR/$casindex_sample/logs/analyse_rare_variants.$logbasename.log"
    
    echo "Running rare variants analysis for sample: $casindex_sample"
    bash "$PIPELINEBASE/common/vcf/wrapper_analyse_rare_variants.sh" \
        -i "$INPUTFILE" \
        -o "$OUTPUTFILE" \
        --pass --mindp 5 --minaltdp 3 --minaltfrac 0.1 --maxsamples 5 \
        --nsssi --freqesp 0.01 --freqexac 0.01 --freqgnomadex 0.01 --freqgnomadge 0.01 \
        --recessive --denovo \
        &>> "$LOGFILE"
done
echo "End: $(date +"%F_%H-%M-%S")"

# Report Generation Step
echo "### Generating VCF Report ###"
echo "Start: $(date +"%F_%H-%M-%S")"
for varcallfile in $(ls "$ANALYSISDIR"/*/*.varcall.list); do
    casindex_sample=$(basename "${varcallfile%%.*}")

    # Clean up previous report files if they exist
    if [ -e "$ANALYSISDIR/$casindex_sample/$casindex_sample.*.report.tsv" ]; then
        rm "$ANALYSISDIR/$casindex_sample/$casindex_sample.*.report.tsv"
    fi

    # Prepare inputs for report generation
    INPUTFILE="$ANALYSISDIR/$casindex_sample/$casindex_sample.rare.vcf"
    OUTPUTFILE="$ANALYSISDIR/$casindex_sample/$casindex_sample.report.tsv"
    LOGFILE="$ANALYSISDIR/$casindex_sample/logs/report_vcf.$logbasename.log"

    echo "Generating report for sample: $casindex_sample"
    bash "$PIPELINEBASE/common/vcf/wrapper_report_vcf.sh" \
        -i "$INPUTFILE" \
        -o "$OUTPUTFILE" \
        &>> "$LOGFILE"
done
echo "End: $(date +"%F_%H-%M-%S")"
