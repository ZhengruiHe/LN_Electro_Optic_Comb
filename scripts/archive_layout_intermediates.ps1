<#
按已审核归档计划移动独立中间稿；不删除文件。默认只检查，-Apply才执行。
每个源和目标都必须位于工作区results/layout内，拒绝符号链接/重解析点。
#>
param(
    [Parameter(Mandatory=$true)][string]$Plan,
    [switch]$Apply
)
$ErrorActionPreference = 'Stop'
$archivePlan = Get-Content -LiteralPath $Plan -Raw | ConvertFrom-Json
$projectPath = (Resolve-Path -LiteralPath $archivePlan.workspace).Path
$layoutPath = [IO.Path]::GetFullPath((Join-Path $projectPath 'results/layout'))
$archivePath = [IO.Path]::GetFullPath($archivePlan.archive)
$prefix = $layoutPath + [IO.Path]::DirectorySeparatorChar
if (-not $archivePath.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) { throw '归档目标越出layout目录' }
$targets = @()
foreach ($entry in $archivePlan.selected) {
    $sourcePath = (Resolve-Path -LiteralPath $entry.source).Path
    $destinationPath = [IO.Path]::GetFullPath($entry.destination)
    if ([IO.Path]::GetDirectoryName($sourcePath) -ne $layoutPath) { throw '源目录不是layout直接子目录' }
    if (-not $destinationPath.StartsWith($archivePath + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) { throw '目标越出指定归档目录' }
    if (Test-Path -LiteralPath $destinationPath) { throw ('目标已存在：' + $destinationPath) }
    $items = @((Get-Item -LiteralPath $sourcePath)) + @(Get-ChildItem -LiteralPath $sourcePath -Recurse -Force)
    if (@($items | Where-Object { $_.Attributes -band [IO.FileAttributes]::ReparsePoint }).Count) { throw '拒绝归档重解析点' }
    foreach ($file in $entry.files) {
        $filePath = [IO.Path]::GetFullPath((Join-Path $sourcePath $file.path))
        if (-not $filePath.StartsWith($sourcePath + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) { throw '文件路径越界' }
        if ((Get-FileHash -LiteralPath $filePath -Algorithm SHA256).Hash.ToLower() -ne $file.sha256) { throw '文件已变化，需要重做计划' }
    }
    $targets += [pscustomobject]@{Source=$sourcePath;Destination=$destinationPath;Entry=$entry}
}
if (-not $Apply) { $targets | Select-Object Source,Destination; exit 0 }
if (-not (Test-Path -LiteralPath $archivePath)) { New-Item -ItemType Directory -Path $archivePath | Out-Null }
foreach ($target in $targets) {
    Move-Item -LiteralPath $target.Source -Destination $target.Destination
    foreach ($file in $target.Entry.files) {
        if ((Get-FileHash -LiteralPath (Join-Path $target.Destination $file.path) -Algorithm SHA256).Hash.ToLower() -ne $file.sha256) { throw '归档后哈希不一致' }
    }
}
[pscustomobject]@{Status='archived_without_deletion';Count=$targets.Count;Archive=$archivePath;AllHashesPreserved=$true} | ConvertTo-Json
