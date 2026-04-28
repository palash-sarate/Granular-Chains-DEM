# sync_drive.ps1
# PowerShell equivalent of Pulse/sync_drive.sh

# Ensure we are in the project root
Set-Location (Split-Path -Parent $PSScriptRoot)

# Configuration
$folders = @("chain_data", "dumping_yard")
$ZIP_DIR = "Zips_dir"

if (!(Test-Path $ZIP_DIR)) {
    New-Item -ItemType Directory -Path $ZIP_DIR -Force | Out-Null
}

Write-Host "Starting optimized sync to Google Drive..." -ForegroundColor Cyan

foreach ($folder in $folders) {
    if (Test-Path $folder) {
        if ($folder -eq "dumping_yard") {
            Write-Host "Using Per-Folder Compression for $folder..." -ForegroundColor Green
            $targetZipDir = Join-Path $ZIP_DIR $folder
            if (!(Test-Path $targetZipDir)) {
                New-Item -ItemType Directory -Path $targetZipDir -Force | Out-Null
            }
            
            # Target leaf simulation directories (depth 2 relative to dumping_yard)
            # equivalent to find "$folder" -mindepth 2 -maxdepth 2 -type d
            $groups = Get-ChildItem -Path $folder -Directory
            foreach ($group in $groups) {
                $sims = Get-ChildItem -Path $group.FullName -Directory
                foreach ($sim in $sims) {
                    # Create a safe filename for the zip
                    $relPath = "$folder/$($group.Name)/$($sim.Name)"
                    $zipName = "$($relPath.Replace('/', '_').Replace('\', '_')).tar.gz"
                    $zipPath = Join-Path $targetZipDir $zipName
                    
                    # Check if compression is needed
                    $needsUpdate = $false
                    if (!(Test-Path $zipPath)) {
                        $needsUpdate = $true
                    } else {
                        # Check if any file in sim is newer than zip
                        $lastWrite = (Get-ChildItem -Path $sim.FullName -Recurse | Measure-Object -Property LastWriteTime -Maximum).Maximum
                        if ($lastWrite -gt (Get-Item $zipPath).LastWriteTime) {
                            $needsUpdate = $true
                        }
                    }
                    
                    if ($needsUpdate) {
                        Write-Host "Compressing $($sim.FullName) -> $zipName" -ForegroundColor Gray
                        tar -czf $zipPath -C (Split-Path $sim.FullName) ($sim.Name)
                    }
                }
            }
            
            Write-Host "Syncing compressed archives to Google Drive..." -ForegroundColor Green
            rclone copy $targetZipDir "gdrive:Granular-Chains-DEM/dumping_yard_zips" `
                --progress `
                --transfers 16 `
                --drive-chunk-size 128M `
                --buffer-size 256M
        }
        else {
            Write-Host "Syncing $folder directly..." -ForegroundColor Green
            rclone copy $folder "gdrive:Granular-Chains-DEM/$folder" `
                --progress `
                --transfers 32 `
                --checkers 64 `
                --fast-list `
                --drive-chunk-size 64M `
                --buffer-size 128M `
                --use-mmap
        }
    }
    else {
        Write-Host "Warning: Folder $folder not found, skipping." -ForegroundColor Yellow
    }
}

Write-Host "Sync complete!" -ForegroundColor Cyan
