#!/bin/bash
# download_ecoli_140.sh
# Download 140 E. coli genomes from NCBI FTP (SynTracker supplementary dataset)
# Usage: bash download_ecoli_140.sh [output_directory]

set -e

OUTPUT_DIR="${1:-./ecoli_140_genomes}"
ACCESSION_FILE="$(dirname "$0")/ecoli_140_accessions.txt"
LOG_FILE="$OUTPUT_DIR/download.log"

mkdir -p "$OUTPUT_DIR"

# Check accession file
if [ ! -f "$ACCESSION_FILE" ]; then
    echo "ERROR: Accession file not found: $ACCESSION_FILE"
    echo "Please ensure ecoli_140_accessions.txt is in the same directory."
    exit 1
fi

# Count total
TOTAL=$(wc -l < "$ACCESSION_FILE" | tr -d ' ')
COMPLETED=0
FAILED=0
SKIPPED=0

echo "========================================"
echo "Downloading $TOTAL E. coli genomes"
echo "Output: $OUTPUT_DIR"
echo "Started at: $(date)"
echo "========================================"

while IFS=$'\t' read -r IDX ACCESSION PHYLOGROUP; do
    OUTPUT="$OUTPUT_DIR/${ACCESSION}.fna.gz"
    
    # Skip if already downloaded and valid
    if [ -f "$OUTPUT" ] && [ -s "$OUTPUT" ]; then
        if gunzip -t "$OUTPUT" 2>/dev/null; then
            SKIPPED=$((SKIPPED + 1))
            continue
        else
            echo "  Removing corrupted file: $OUTPUT"
            rm -f "$OUTPUT"
        fi
    fi
    
    # Parse accession to build URL
    # GCA_001575755.1 -> prefix=GCA, number=001575755, version=1
    BASE=${ACCESSION%%.*}
    PREFIX=$(echo "$BASE" | cut -d'_' -f1)
    NUMBER=$(echo "$BASE" | cut -d'_' -f2)
    
    SEG1=${NUMBER:0:3}
    SEG2=${NUMBER:3:3}
    SEG3=${NUMBER:6:3}
    
    # Build assembly directory name
    # GCA_001575755.1 -> GCA_001575755.1_ASM157575v1
    # Remove leading zeros from number for ASM part
    ASM_NUM=$(echo "$NUMBER" | sed 's/^0*//')
    [ -z "$ASM_NUM" ] && ASM_NUM="0"
    
    ASM_DIR="${ACCESSION}_ASM${ASM_NUM}v1"
    URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${PREFIX}/${SEG1}/${SEG2}/${SEG3}/${ASM_DIR}/${ASM_DIR}_genomic.fna.gz"
    
    echo "[$((COMPLETED + FAILED + SKIPPED + 1))/$TOTAL] Downloading $ACCESSION ($PHYLOGROUP)..."
    
    SUCCESS=0
    for ATTEMPT in 1 2; do
        if [ "$ATTEMPT" = "2" ]; then
            # Try v2 version
            ASM_DIR="${ACCESSION}_ASM${ASM_NUM}v2"
            URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${PREFIX}/${SEG1}/${SEG2}/${SEG3}/${ASM_DIR}/${ASM_DIR}_genomic.fna.gz"
            echo "  Retrying with v2: $ASM_DIR"
        fi
        
        # Download with resume support and 5-minute timeout per attempt
        if curl -s -L -C - --max-time 300 -o "$OUTPUT" "$URL" 2>/dev/null; then
            # Verify it's a valid gzip/FASTA
            if [ -f "$OUTPUT" ] && [ -s "$OUTPUT" ]; then
                if gunzip -t "$OUTPUT" 2>/dev/null; then
                    # Check first line is FASTA header
                    FIRST=$(gunzip -c "$OUTPUT" 2>/dev/null | head -1)
                    if echo "$FIRST" | grep -q '^>'; then
                        SUCCESS=1
                        COMPLETED=$((COMPLETED + 1))
                        break
                    else
                        echo "  Warning: File exists but not valid FASTA, retrying..."
                        rm -f "$OUTPUT"
                    fi
                else
                    echo "  Warning: Invalid gzip, retrying..."
                    rm -f "$OUTPUT"
                fi
            fi
        fi
    done
    
    if [ "$SUCCESS" -eq 0 ]; then
        echo "[FAIL] $ACCESSION" | tee -a "$LOG_FILE"
        FAILED=$((FAILED + 1))
        rm -f "$OUTPUT"
    else
        echo "[OK] $ACCESSION" >> "$LOG_FILE"
    fi
    
    # Progress every 10 genomes
    if [ $(( (COMPLETED + FAILED + SKIPPED) % 10 )) -eq 0 ] && [ "$COMPLETED" -gt 0 ]; then
        echo "  Progress: $COMPLETED completed, $FAILED failed, $SKIPPED skipped"
    fi
done < "$ACCESSION_FILE"

echo ""
echo "========================================"
echo "Download Complete"
echo "Completed: $COMPLETED"
echo "Failed: $FAILED"
echo "Skipped (already exist): $SKIPPED"
echo "Total: $TOTAL"
echo "Finished at: $(date)"
echo "========================================"

if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "Failed downloads logged in: $LOG_FILE"
    echo "You can retry failed downloads by running this script again."
fi
