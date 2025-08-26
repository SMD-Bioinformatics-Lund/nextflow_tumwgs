#!/bin/bash

# Function to write log messages
function log_message {
    local message=$1
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $message" >> "$LOG_FILE"
}

# Function to print error messages, log them, and exit
function error_exit {
    echo "Error: $1" >&2
    log_message "ERROR: $1"
    exit 1
}

# Function to copy a file
function copy_file {
    local source_file="$1"
    local dest_file="$2"
    scp "$source_file" "$DEST_HOST:$dest_file"
    if [ $? -eq 0 ]; then
        log_message "Successfully copied $source_file to $DEST_HOST:$dest_file"
    else
        error_exit "Failed to copy $source_file to $DEST_HOST:$dest_file"
    fi
}

# Function to copy a directory
function copy_directory {
    local source_dir="$1"
    local dest_dir="$2"
    scp -r "$source_dir" "$DEST_HOST:$dest_dir"
    if [ $? -eq 0 ]; then
        log_message "Successfully copied directory $source_dir to $DEST_HOST:$dest_dir"
    else
        error_exit "Failed to copy directory $source_dir to $DEST_HOST:$dest_dir"
    fi
}

# Get the current script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" || error_exit "Failed to determine script directory"
LOG_FILE="$DIR/transfer.log"
echo "$DIR"

# Define destination host and pipeline path 
DEST_HOST="rs-fs1.lunarc.lu.se"
PIPELINE_DEST="/fs1/saile/prj/pipeline_test/nextflow_tumwgs"

# Ensure the log file exists and is writable
touch "$LOG_FILE" || error_exit "Cannot create or write to log file $LOG_FILE"

# Main copy operations
log_message "Starting file transfer operations."

echo "Copying the main.nf files..."
copy_file "$DIR/main.nf" "$PIPELINE_DEST"

echo "Copying configuration file..."
copy_file "$DIR/configs/nextflow.hopper.config" "$PIPELINE_DEST/nextflow.config"
copy_directory "$DIR/configs" "$PIPELINE_DEST"

echo "Copying other files..."
copy_directory "$DIR/resources" "$PIPELINE_DEST"
copy_directory "$DIR/bin" "$PIPELINE_DEST"
copy_directory "$DIR/modules" "$PIPELINE_DEST"
copy_directory "$DIR/subworkflows" "$PIPELINE_DEST"
copy_directory "$DIR/workflows" "$PIPELINE_DEST"

log_message "All files copied successfully."
echo "All files copied successfully."