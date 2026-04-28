#!/bin/bash

# Ensure we are in the project root
cd "$(dirname "$0")/.." || exit 1

# Configuration
folders=("chain_data" "dumping_yard")
ZIP_DIR="Zips_dir"
mkdir -p "$ZIP_DIR"

echo "Starting download from Google Drive..."

for folder in "${folders[@]}"; do
    if [[ "$folder" == "dumping_yard" ]]; then
        echo "Downloading compressed archives for $folder..."
        mkdir -p "$ZIP_DIR/$folder"
        
        # rclone copy only adds/updates files (does NOT delete anything locally)
        rclone copy "gdrive:Granular-Chains-DEM/dumping_yard_zips" "$ZIP_DIR/$folder" \
            --progress \
            --transfers 16 \
            --drive-chunk-size 128M \
            --buffer-size 256M

        echo "Extracting compressed archives..."
        for zip_file in "$ZIP_DIR/$folder"/*.tar.gz; do
            if [ -f "$zip_file" ]; then
                # Extract files. tar will overwrite if files are newer or different.
                # No deletion of local files that aren't in the zip.
                echo "Extracting $zip_file"
                tar -xzf "$zip_file"
            fi
        done
    else
        echo "Downloading $folder directly..."
        # rclone copy only adds/updates files (does NOT delete anything locally)
        rclone copy "gdrive:Granular-Chains-DEM/$folder" "$folder" \
            --progress \
            --transfers 32 \
            --checkers 64 \
            --fast-list \
            --drive-chunk-size 64M \
            --buffer-size 128M \
            --use-mmap
    fi
done

echo "Download complete!"
