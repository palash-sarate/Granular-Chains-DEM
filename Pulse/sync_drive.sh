#!/bin/bash
# Sync local project folders to Google Drive

# 1. Ensure we are in the project root
PROJECT_ROOT="/home/guest/palash/Granular-Chains-DEM"
cd $PROJECT_ROOT

# 2. Check if rclone is configured
if ! rclone listremotes | grep -q "gdrive:"; then
    echo "Error: 'gdrive' remote not found. Please run 'rclone config' first."
    exit 1
fi

echo "Starting sync to Google Drive..."

# Sync folders
# We use 'copy' or 'sync'. 'sync' makes the destination match the source (deletes files on Drive if deleted locally).
# 'copy' only adds/updates files. I'll use 'copy' for safety unless you prefer 'sync'.

folders=("chain_data" "dumping_yard" "PBS_Output")

for folder in "${folders[@]}"; do
    if [ -d "$folder" ]; then
        echo "Syncing $folder..."
        rclone copy "$folder" "gdrive:Granular-Chains-DEM/$folder" --progress
    else
        echo "Warning: Folder $folder not found, skipping."
    fi
done

echo "Sync complete!"
