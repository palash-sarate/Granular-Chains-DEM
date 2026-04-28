# download_drive.ps1
# PowerShell equivalent of Pulse/download_drive.sh

# Ensure we are in the project root
Set-Location (Split-Path -Parent $PSScriptRoot)

# Configuration
$folders = @("chain_data", "dumping_yard")
$ZIP_DIR = "Zips_dir"

if (!(Test-Path $ZIP_DIR)) {
    New-Item -ItemType Directory -Path $ZIP_DIR -Force | Out-Null
}

Write-Host "Starting download from Google Drive..." -ForegroundColor Cyan

foreach ($folder in $folders) {
    if ($folder -eq "dumping_yard") {
        Write-Host "Downloading compressed archives for $folder..." -ForegroundColor Green
        $targetZipDir = Join-Path $ZIP_DIR $folder
        if (!(Test-Path $targetZipDir)) {
            New-Item -ItemType Directory -Path $targetZipDir -Force | Out-Null
        }
        
        # rclone copy only adds/updates files
        rclone copy "gdrive:Granular-Chains-DEM/dumping_yard_zips" $targetZipDir `
            --progress `
            --transfers 16 `
            --drive-chunk-size 128M `
            --buffer-size 256M

        Write-Host "Extracting compressed archives..." -ForegroundColor Cyan
        $zips = Get-ChildItem -Path $targetZipDir -Filter "*.tar.gz"
        foreach ($zip in $zips) {
            Write-Host "Extracting $($zip.Name)" -ForegroundColor Gray
            # Windows 10/11 has tar built-in
            tar -xzf $zip.FullName
        }
    }
    else {
        Write-Host "Downloading $folder directly..." -ForegroundColor Green
        # rclone copy only adds/updates files
        rclone copy "gdrive:Granular-Chains-DEM/$folder" $folder `
            --progress `
            --transfers 32 `
            --checkers 64 `
            --fast-list `
            --drive-chunk-size 64M `
            --buffer-size 128M `
            --use-mmap
    }
}

Write-Host "Download complete!" -ForegroundColor Cyan
