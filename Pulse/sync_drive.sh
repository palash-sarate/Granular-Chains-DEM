#!/bin/bash

# Ensure we are in the project root
cd "$(dirname "$0")/.." || exit 1

# Configuration
folders=("chain_data" "dumping_yard")
ZIP_DIR="Zips_dir"
mkdir -p "$ZIP_DIR"

echo "Starting optimized sync to Google Drive..."

for folder in "${folders[@]}"; do
    if [ -d "$folder" ]; then
        if [[ "$folder" == "dumping_yard" ]]; then
            echo "Using Per-Folder Compression for $folder..."
            mkdir -p "$ZIP_DIR/$folder"
            
            # Target leaf simulation directories (depth 2 relative to dumping_yard)
            # e.g., dumping_yard/Relax_3d_Library_Gen/relax_N4_state_0
            find "$folder" -mindepth 2 -maxdepth 2 -type d | while read -r sub; do
                # Create a safe filename for the zip (replace / with _)
                ZIP_NAME=$(echo "$sub" | tr '/' '_').tar.gz
                ZIP_PATH="$ZIP_DIR/$folder/$ZIP_NAME"
                
                # Compress if zip doesn't exist OR if any file in folder is newer than zip
                if [ ! -f "$ZIP_PATH" ] || [ "$(find "$sub" -newer "$ZIP_PATH" | wc -l)" -gt 0 ]; then
                    echo "Compressing $sub -> $ZIP_NAME"
                    tar -czf "$ZIP_PATH" "$sub"
                fi
            done
            
            echo "Syncing compressed archives to Google Drive..."
            # rclone copy only adds/updates files (does NOT delete anything on remote)
            rclone copy "$ZIP_DIR/$folder" "gdrive:Granular-Chains-DEM/dumping_yard_zips" \
                --progress \
                --transfers 16 \
                --drive-chunk-size 128M \
                --buffer-size 256M
        else
            # Standard high-parallelism sync for other folders (like chain_data)
            echo "Syncing $folder directly..."
            # rclone copy only adds/updates files (does NOT delete anything on remote)
            rclone copy "$folder" "gdrive:Granular-Chains-DEM/$folder" \
                --progress \
                --transfers 32 \
                --checkers 64 \
                --fast-list \
                --drive-chunk-size 64M \
                --buffer-size 128M \
                --use-mmap
        fi
    else
        echo "Warning: Folder $folder not found, skipping."
    fi
done

echo "Sync complete!"
