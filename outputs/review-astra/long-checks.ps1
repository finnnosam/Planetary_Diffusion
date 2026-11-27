# Run one requested diagnostic at a time. Each invocation creates a fresh output directory.
param([Parameter(Mandatory=$true)][ValidateSet('Repeat','Chain','LongSeed2')][string]$Check)
$ErrorActionPreference = 'Stop'
$root = 'E:\projects\Planetary_Diffusion'
Set-Location -LiteralPath $root
$py = Join-Path $root '.venv\Scripts\python.exe'
$env:NUMBA_NUM_THREADS = '2'
$destination = Join-Path $root ("outputs\review-astra\runs\" + $Check + "-" + [guid]::NewGuid().ToString('N').Substring(0,8))
New-Item -ItemType Directory -Path $destination | Out-Null
Write-Host "Results in $destination"
$plates = Join-Path $destination 'plates'
$common = @('-u','-m','scripts.simulate_plates','--n','217','--out',$plates,'--param','threads=2','--param','rift_traction=True')
if ($Check -eq 'LongSeed2') {
    & $py @common --seeds 1002 --years 500 --save-every 10
    exit $LASTEXITCODE
}
if ($Check -eq 'Repeat') {
    & $py @common --seeds 1001 --years 250 --save-every 5
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    $original = Join-Path $root 'outputs\planet-calibration\sweep-20261122-150810\combo-01\seed1001'
    $repeat = Join-Path $plates 'seed1001'
    & $py (Join-Path $root 'tools\compare_runs.py') $original $repeat |
        Tee-Object -FilePath (Join-Path $destination 'shared-values.txt')
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    $log = Get-Content -Raw (Join-Path $destination 'shared-values.txt')
    if ($log -notmatch '(?m)^IDENTICAL ' -or $log -match '(?m)^  only in [AB]:') { exit 1 }
    # The repository comparator reads array values and stats, excluding timing and CLI overrides.
    # Check snapshot inventory and the plate-rotation tables as well.
    $left = Join-Path $original 'snapshots'
    $right = Join-Path $repeat 'snapshots'
    $namesA = @(Get-ChildItem -LiteralPath $left -File -Name | Sort-Object)
    $namesB = @(Get-ChildItem -LiteralPath $right -File -Name | Sort-Object)
    if (Compare-Object -ReferenceObject $namesA -DifferenceObject $namesB) { throw 'Snapshot inventories differ' }
    foreach ($name in ($namesA | Where-Object { $_ -like '*-plates.json' })) {
        if ((Get-FileHash -LiteralPath (Join-Path $left $name)).Hash -ne (Get-FileHash -LiteralPath (Join-Path $right $name)).Hash) {
            throw "Plate table differs: $name"
        }
    }
    Write-Host 'Snapshot fields, plate tables and recorded scientific statistics agree.'
    exit 0
}
& $py @common --seeds 1001 --years 250 --save-last 40 --save-every 1
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
$landscape = Join-Path $destination 'landscape'
& $py -u -m scripts.run_planet_landscape --plates (Join-Path $plates 'seed1001') --out $landscape
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $py -u -m scripts.planet_meso_check --landscape $landscape --out (Join-Path $destination 'meso')
exit $LASTEXITCODE
