#!/bin/bash

# Base directory to search for JSON files
BASE_DIR="/project/3018051.01/ruggero/bids/"

# Find all .json files in subdirectories
find "$BASE_DIR" -type f -name "*.json" | while read -r file; do
    # Process the B0FieldIdentifier field
    jq 'if has("B0FieldIdentifier") and (.B0FieldIdentifier | type == "string") then .B0FieldIdentifier |= gsub("<<|>>"; "") else . end' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
    
    # Process the B0FieldSource field
    jq 'if has("B0FieldSource") and (.B0FieldSource | type == "string") then .B0FieldSource |= gsub("<<|>>"; "") else . end' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"

    # Check for errors
    if [[ $? -eq 0 ]]; then
        echo "Processed: $file"
    else
        echo "Error processing: $file"
    fi
done

echo "Processing complete."

