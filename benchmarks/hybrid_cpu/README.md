# Hybrid-CPU last-bit reproducer (Unit 41)

Investigation only: no fix and no change to any result (DECISIONS.md, Unit 41).

`repro.cpp` uses Eigen alone. Build it from a configured release tree, whose vcpkg install holds
Eigen, in a Visual Studio x64 developer prompt, with the project's floating-point flags:

```text
cl /nologo /O2 /DNDEBUG /EHsc /std:c++latest /utf-8 /permissive- /fp:precise /MD ^
   /I build\release\vcpkg_installed\x64-windows-static-md\include\eigen3 ^
   benchmarks\hybrid_cpu\repro.cpp /Fe:repro.exe
```

Run it pinned to one logical CPU at a time (`/affinity` takes a hex mask: 1 is CPU 0, 4 is CPU 2):

```text
start /b /wait /affinity 1 repro.exe auto
start /b /wait /affinity 4 repro.exe auto
start /b /wait /affinity 4 repro.exe fixed 49152 3145728 20971520
start /b /wait /affinity 1 repro.exe fixed 32768 3145728 20971520
```

On an Intel Core Ultra 5 225F (CPUs 0, 1, 6-9 P cores: L1d 48 KiB, L2 3 MiB; CPUs 2-5 E cores:
L1d 32 KiB, L2 4 MiB), `auto` prints one SparseLU hash on the P cores and another on the E cores.
With the cache sizes fixed, every core prints the same hash. Fixing only the L1 size moves a P
core to the E-core result and back; the L2 size has no effect.

`diode_2d.deck` shows the same thing through `nitcad run`. Its current table differs in the last
bits between a run pinned to a P core and one pinned to an E core.
