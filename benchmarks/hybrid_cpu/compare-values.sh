#!/bin/sh
# compare-values.sh <dir A> <dir B>: per traced entry, compares two value traces (the Unit 41
# measurement listener) block by block (test case and section). A block whose sequence of calls
# (tag and length per line) differs ran a different number of evaluations: iterations, steps or
# cuts changed. A block with the same sequence but different hashes changed values only; for the
# result.* tags the size is reported as max |a - b| over max |a| of each array.
A="$1"; B="$2"
for t in "$A"/*.trace; do
  n=$(basename "$t")
  awk -v F2="$B/$n" -v NAME="$n" '
    function flush(   i, a, b, na, nb, same, k, va, vb, m, d, j, x, y, tag) {
      if (n1 == 0 && n2 == 0) return
      blocks++
      same = (n1 == n2)
      for (i = 1; same && i <= n1; i++) {
        split(b1[i], a, "\t"); split(b2[i], b, "\t")
        if (a[2] != b[2] || a[3] != b[3]) same = 0
      }
      if (!same) {
        seq++
        if (seq <= 6) seqlist = seqlist "\n    " head " (" n1 " against " n2 " calls)"
      } else {
        changed = 0
        for (i = 1; i <= n1; i++) {
          split(b1[i], a, "\t"); split(b2[i], b, "\t")
          if (a[4] == b[4]) continue
          changed = 1
          tag = a[2]
          if (tag !~ /^result\./) { internal[tag]++; continue }
          na = split(a[5], va, " "); nb = split(b[5], vb, " ")
          m = 0; d = 0
          for (j = 1; j <= na; j++) {
            x = va[j] + 0; y = vb[j] + 0
            if (x < 0) x = -x
            if (x > m) m = x
            x = va[j] - vb[j]; if (x < 0) x = -x
            if (x > d) d = x
          }
          r = (m > 0) ? d / m : d
          if (!(tag in worst) || r > worst[tag]) worst[tag] = r
        }
        if (changed) valued++
      }
      n1 = 0; n2 = 0
    }
    function read2(   l) {
      n2 = 0
      while ((getline l < F2) > 0) {
        if (l ~ /^(case|section)\t/) { next2 = l; return }
        b2[++n2] = l
      }
      next2 = ""
    }
    BEGIN { n1 = 0; n2 = 0; blocks = 0; seq = 0; valued = 0; started = 0 }
    /^(case|section)\t/ {
      if (started) { read2(); flush() } else { while ((getline l < F2) > 0) if (l ~ /^(case|section)\t/) break }
      started = 1
      if ($0 ~ /^case\t/) casename = substr($0, 6)
      head = casename (($0 ~ /^section\t/) ? " / " substr($0, 9) : "")
      next
    }
    { b1[++n1] = $0 }
    END {
      read2(); flush()
      printf "%s: %d blocks, %d with a different call sequence, %d with different values only\n", NAME, blocks, seq, valued
      if (seq > 0) printf "  sequence changes:%s\n", seqlist
      for (k in worst) printf "  %-36s largest relative difference %.3g\n", k, worst[k]
      for (k in internal) ni++
      if (ni > 0) printf "  (%d assembler tags differ as well)\n", ni
    }' "$t"
done
