# Makes the NiTCAD icon from the source art: the wafer alone on a transparent background.
# Usage: powershell -File make_icon.ps1 -Src nitcad_source.png -Png nitcad.png -Ico nitcad.ico
# The wafer's top face is fitted as an ellipse (moments of its strongly coloured pixels, holes
# filled), then grown by the rim ring and extended down by the dark side wall. Everything outside
# that shape - the white tile, the drop shadow and the wordmark - becomes transparent; the edge is
# anti-aliased by 4x4 supersampling.
param([string]$Src, [string]$Png, [string]$Ico)
Add-Type -AssemblyName System.Drawing
Add-Type -ReferencedAssemblies System.Drawing -TypeDefinition @'
using System; using System.Collections.Generic; using System.Drawing; using System.Drawing.Imaging;
using System.Runtime.InteropServices;
public static class WaferIcon {
    public static string Report = "";
    public static Bitmap Cutout(Bitmap src, int textTop, int side, double ring) {
        int w = src.Width, h = src.Height;
        var d = src.LockBits(new Rectangle(0, 0, w, h), ImageLockMode.ReadOnly, PixelFormat.Format32bppArgb);
        var px = new byte[w * h * 4]; Marshal.Copy(d.Scan0, px, 0, px.Length); src.UnlockBits(d);
        int fh = Math.Min(h, textTop);
        // Face pixels: strongly coloured (cells and grid lines); the shadow and the rim are not.
        var face = new bool[w * fh];
        for (int i = 0; i < w * fh; i++) {
            int b = px[4 * i], g = px[4 * i + 1], r = px[4 * i + 2];
            int mx = Math.Max(r, Math.Max(g, b)), mn = Math.Min(r, Math.Min(g, b));
            face[i] = px[4 * i + 3] > 200 && mx - mn >= 110;
        }
        // Fill the holes (pale cells): whatever the outside can't reach is face.
        var outside = new bool[w * fh]; var q = new Queue<int>();
        for (int x = 0; x < w; x++) { q.Enqueue(x); q.Enqueue((fh - 1) * w + x); }
        for (int y = 0; y < fh; y++) { q.Enqueue(y * w); q.Enqueue(y * w + w - 1); }
        while (q.Count > 0) {
            int i = q.Dequeue(); if (outside[i] || face[i]) continue; outside[i] = true;
            int x = i % w, y = i / w;
            if (x > 0) q.Enqueue(i - 1); if (x < w - 1) q.Enqueue(i + 1);
            if (y > 0) q.Enqueue(i - w); if (y < fh - 1) q.Enqueue(i + w);
        }
        double n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        for (int i = 0; i < w * fh; i++) if (!outside[i]) { double x = i % w, y = i / w; n++; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; }
        double cx = sx / n, cy = sy / n, vxx = sxx / n - cx * cx, vyy = syy / n - cy * cy, vxy = sxy / n - cx * cy;
        double tr = vxx + vyy, det = vxx * vyy - vxy * vxy, disc = Math.Sqrt(tr * tr / 4 - det);
        double th = 0.5 * Math.Atan2(2 * vxy, vxx - vyy);
        double a = 2 * Math.Sqrt(tr / 2 + disc) + ring, bb = 2 * Math.Sqrt(tr / 2 - disc) + ring;
        double ct = Math.Cos(th), st = Math.Sin(th);
        Func<double, double, double> E = (x, y) => {
            double dx = x - cx, dy = y - cy, u = dx * ct + dy * st, v = -dx * st + dy * ct;
            return u * u / (a * a) + v * v / (bb * bb);
        };
        // Side wall: below the face's lowest point at the centre column, the dark band before the
        // shadow lightens.
        int xc = (int)Math.Round(cx), ybot = (int)cy;
        while (ybot < fh - 1 && E(xc, ybot + 1) <= 1) ybot++;
        int ywall = ybot;
        for (int y = ybot; y < Math.Min(fh, ybot + 80); y++) {
            int i = 4 * (y * w + xc); int mx = Math.Max(px[i], Math.Max(px[i + 1], px[i + 2]));
            if (mx < 60) ywall = y;
        }
        double t = ywall - ybot + 1.0;
        Report = string.Format("face centre ({0:F1},{1:F1}) semi-axes {2:F1} x {3:F1} angle {4:F2} deg, wall {5:F0} px",
                               cx, cy, a, bb, th * 180 / Math.PI, t);
        // Inside: the face ellipse swept straight down by 0..t (the wall).
        Func<double, double, bool> inside = (x, y) => {
            for (int k = 0; k <= 32; k++) if (E(x, y - t * k / 32) <= 1) return true;
            return false;
        };
        int x0 = w, y0 = fh, x1 = 0, y1 = 0;
        var o = new byte[w * fh * 4];
        for (int y = 0; y < fh; y++) for (int x = 0; x < w; x++) {
            // Far from the edge, skip the supersampling.
            double e = E(x + 0.5, y + 0.5);
            int cov = 0;
            if (e > 1.6 && E(x + 0.5, y + 0.5 - t) > 1.6) cov = 0;
            else
                for (int s2 = 0; s2 < 16; s2++)
                    if (inside(x + (s2 % 4 + 0.5) / 4, y + (s2 / 4 + 0.5) / 4)) cov++;
            if (cov == 0) continue;
            int i = 4 * (y * w + x);
            o[i] = px[i]; o[i + 1] = px[i + 1]; o[i + 2] = px[i + 2]; o[i + 3] = (byte)(255 * cov / 16);
            x0 = Math.Min(x0, x); y0 = Math.Min(y0, y); x1 = Math.Max(x1, x); y1 = Math.Max(y1, y);
        }
        var cut = new Bitmap(w, fh, PixelFormat.Format32bppArgb);
        var cd = cut.LockBits(new Rectangle(0, 0, w, fh), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
        Marshal.Copy(o, 0, cd.Scan0, o.Length); cut.UnlockBits(cd);
        // Centre the wafer on a transparent square, spanning 96 % of the width.
        var res = new Bitmap(side, side, PixelFormat.Format32bppArgb);
        using (var gr = Graphics.FromImage(res)) {
            gr.Clear(Color.Transparent);
            gr.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
            gr.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.HighQuality;
            int bw = x1 - x0 + 1, bh = y1 - y0 + 1;
            float s = side * 0.96f / Math.Max(bw, bh), dw = bw * s, dh = bh * s;
            gr.DrawImage(cut, new RectangleF((side - dw) / 2, (side - dh) / 2, dw, dh),
                         new RectangleF(x0, y0, bw, bh), GraphicsUnit.Pixel);
        }
        cut.Dispose(); return res;
    }
    // Halving steps keep every source pixel in the average down to small sizes.
    public static Bitmap Shrink(Bitmap big, int s) {
        Bitmap cur = big;
        while (cur.Width > s) {
            int n = Math.Max(s, cur.Width / 2);
            var nb = new Bitmap(n, n, PixelFormat.Format32bppArgb);
            using (var g = Graphics.FromImage(nb)) {
                g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
                g.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.HighQuality;
                g.CompositingMode = System.Drawing.Drawing2D.CompositingMode.SourceCopy;
                g.DrawImage(cur, 0, 0, n, n);
            }
            if (cur != big) cur.Dispose(); cur = nb;
        }
        return cur == big ? (Bitmap)big.Clone() : cur;
    }
}
'@
$img = New-Object System.Drawing.Bitmap $Src
# Rows from 740 down hold the wordmark; the rim ring outside the coloured face is about 4 px.
$master = [WaferIcon]::Cutout($img, 740, 1024, 4.0)
$img.Dispose()
[WaferIcon]::Report
$master.Save($Png, [System.Drawing.Imaging.ImageFormat]::Png)
$sizes = 16, 20, 24, 32, 40, 48, 64, 128, 256
$blobs = foreach ($s in $sizes) {
    $b = [WaferIcon]::Shrink($master, $s); $ms = New-Object System.IO.MemoryStream
    $b.Save($ms, [System.Drawing.Imaging.ImageFormat]::Png); $b.Dispose(); , $ms.ToArray()
}
$master.Dispose()
$out = New-Object System.IO.MemoryStream; $w = New-Object System.IO.BinaryWriter $out
$w.Write([UInt16]0); $w.Write([UInt16]1); $w.Write([UInt16]$sizes.Count)
$offset = 6 + 16 * $sizes.Count
for ($k = 0; $k -lt $sizes.Count; $k++) {
    $d = if ($sizes[$k] -ge 256) { 0 } else { $sizes[$k] }
    $w.Write([byte]$d); $w.Write([byte]$d); $w.Write([byte]0); $w.Write([byte]0)
    $w.Write([UInt16]1); $w.Write([UInt16]32); $w.Write([UInt32]$blobs[$k].Length); $w.Write([UInt32]$offset)
    $offset += $blobs[$k].Length
}
foreach ($bl in $blobs) { $w.Write([byte[]]$bl) }
$w.Flush(); [System.IO.File]::WriteAllBytes($Ico, $out.ToArray())
"icon: $($sizes -join ',') px"
