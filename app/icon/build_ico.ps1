# Builds nitcad.ico: the small wafer-only PNGs plus the original ICO's 64/128/256 entries.
# Usage: powershell -File build_ico.ps1 -OrigIco nitcad.ico -SmallDir <scratch dir> -Out nitcad.ico
# (keeps the 64/128/256 entries of OrigIco, which show the full design with the wordmark)
param([string]$OrigIco, [string]$SmallDir, [string]$Out)
$orig = [System.IO.File]::ReadAllBytes($OrigIco)
$entries = New-Object System.Collections.Generic.List[object]
foreach ($s in 16, 20, 24, 32, 40, 48) {
    $entries.Add(@{ Size = $s; Data = [System.IO.File]::ReadAllBytes((Join-Path $SmallDir "small_$s.png")) })
}
$count = [BitConverter]::ToUInt16($orig, 4)
for ($k = 0; $k -lt $count; $k++) {
    $e = 6 + 16 * $k
    $size = if ($orig[$e] -eq 0) { 256 } else { $orig[$e] }
    if ($size -lt 64) { continue }
    $len = [BitConverter]::ToInt32($orig, $e + 8); $off = [BitConverter]::ToInt32($orig, $e + 12)
    $data = New-Object byte[] $len; [Array]::Copy($orig, $off, $data, 0, $len)
    $entries.Add(@{ Size = $size; Data = $data })
}
$ms = New-Object System.IO.MemoryStream; $w = New-Object System.IO.BinaryWriter $ms
$w.Write([UInt16]0); $w.Write([UInt16]1); $w.Write([UInt16]$entries.Count)
$offset = 6 + 16 * $entries.Count
foreach ($en in $entries) {
    $b = if ($en.Size -ge 256) { 0 } else { $en.Size }
    $w.Write([byte]$b); $w.Write([byte]$b); $w.Write([byte]0); $w.Write([byte]0)
    $w.Write([UInt16]1); $w.Write([UInt16]32); $w.Write([UInt32]$en.Data.Length); $w.Write([UInt32]$offset)
    $offset += $en.Data.Length
}
foreach ($en in $entries) { $w.Write([byte[]]$en.Data) }
$w.Flush(); [System.IO.File]::WriteAllBytes($Out, $ms.ToArray())
"wrote $($entries.Count) entries: $(($entries | ForEach-Object { $_.Size }) -join ',')"
