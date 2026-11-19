# Renders the small NiTCAD icon sizes: the wafer alone on the light tile (no wordmark).
# Usage: powershell -File render_small.ps1 -Src nitcad.png -OutDir <scratch dir>
param([string]$Src, [string]$OutDir)
Add-Type -AssemblyName System.Drawing
$img = New-Object System.Drawing.Bitmap $Src
$crop = New-Object System.Drawing.Rectangle 100, 100, 852, 636
function Get-Rounded([single]$x, [single]$y, [single]$w, [single]$h, [single]$r) {
    $p = New-Object System.Drawing.Drawing2D.GraphicsPath
    $d = 2 * $r
    $p.AddArc($x, $y, $d, $d, 180, 90); $p.AddArc($x + $w - $d, $y, $d, $d, 270, 90)
    $p.AddArc($x + $w - $d, $y + $h - $d, $d, $d, 0, 90); $p.AddArc($x, $y + $h - $d, $d, $d, 90, 90)
    $p.CloseFigure(); return $p
}
foreach ($s in 16, 20, 24, 32, 40, 48) {
    # Draw at 8x, then shrink in halving steps so every source pixel contributes.
    $big = 8 * $s
    $bmp = New-Object System.Drawing.Bitmap $big, $big, ([System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.SmoothingMode = 'AntiAlias'; $g.InterpolationMode = 'HighQualityBicubic'; $g.PixelOffsetMode = 'HighQuality'
    $g.Clear([System.Drawing.Color]::Transparent)
    $tile = Get-Rounded 0 0 ($big - 1) ($big - 1) ($big * 0.18)
    $g.FillPath((New-Object System.Drawing.SolidBrush ([System.Drawing.Color]::FromArgb(255, 249, 251, 253))), $tile)
    $g.DrawPath((New-Object System.Drawing.Pen ([System.Drawing.Color]::FromArgb(255, 200, 208, 220)), ($big / 16)), $tile)
    $inset = $big * 0.07
    $w = $big - 2 * $inset; $h = $w * $crop.Height / $crop.Width
    $dest = New-Object System.Drawing.RectangleF $inset, (($big - $h) / 2), $w, $h
    $g.DrawImage($img, $dest, $crop, [System.Drawing.GraphicsUnit]::Pixel)
    $g.Dispose()
    $cur = $bmp
    while ($cur.Width -gt $s) {
        $n = [Math]::Max($s, [int]($cur.Width / 2))
        $nb = New-Object System.Drawing.Bitmap $n, $n, ([System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
        $g = [System.Drawing.Graphics]::FromImage($nb)
        $g.InterpolationMode = 'HighQualityBicubic'; $g.PixelOffsetMode = 'HighQuality'; $g.CompositingMode = 'SourceCopy'
        $g.DrawImage($cur, 0, 0, $n, $n); $g.Dispose(); $cur.Dispose(); $cur = $nb
    }
    $cur.Save((Join-Path $OutDir "small_$s.png"), [System.Drawing.Imaging.ImageFormat]::Png); $cur.Dispose()
}
$img.Dispose()
