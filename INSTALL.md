# Installing NiTCAD on a new machine

NiTCAD is built from source. It is Windows-only (x64) and needs the MSVC toolchain, so there is no Linux or macOS
build. These steps take a clean Windows 10/11 machine to a working `nitcad.exe` and `nitcad-app.exe`.

Allow about 10 GB of disk and some time for the first build: vcpkg compiles Eigen and Catch2 once, then caches them.

## 1. Install the toolchain

1. **Git**: <https://git-scm.com/download/win>, or `winget install Git.Git`.
2. **Visual Studio 2026** (Community is enough) or the matching **Build Tools**. In the installer select the
   *Desktop development with C++* workload, then under *Individual components* make sure these are ticked:
   - MSVC toolset **14.51** (x64/x86 build tools)
   - Windows SDK **10.0.26100.0**
3. **CMake 4.3 or later**: `winget install Kitware.CMake`
4. **Ninja**: `winget install Ninja-build.Ninja`

Close and reopen your terminal afterwards so the new `PATH` entries are picked up.

## 2. Install the pinned vcpkg

The build refuses any vcpkg checkout that is not at the commit named by `builtin-baseline` in `vcpkg.json`
(`fbb0f7bb200b07a9eb9081c7a3cf51d1aa1c51a1`). Do not use the vcpkg bundled with Visual Studio.

```bat
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
git -C C:\vcpkg checkout fbb0f7bb200b07a9eb9081c7a3cf51d1aa1c51a1
C:\vcpkg\bootstrap-vcpkg.bat -disableMetrics
```

## 3. Get the source

```bat
git clone https://github.com/nihaltodurkar-cms/NiTCAD.git
cd NiTCAD
```

## 4. Build and test

Open **x64 Native Tools Command Prompt for VS** (or *Developer Command Prompt* with the x64 environment), then:

```bat
rem The developer prompt points VCPKG_ROOT at Visual Studio's own vcpkg; override it.
set VCPKG_ROOT=C:\vcpkg

cmake --preset release
cmake --build --preset release
ctest --preset release
```

The first configure downloads and builds Eigen 5.0.1 and Catch2 3.16.0 through vcpkg. Output goes to
`build\release\`. Use `debug` instead of `release` for a Debug build (slower, full test run takes longer).
Warnings are errors (`/W4 /WX`).

## 5. Check it works

```bat
build\release\nitcad.exe run examples\pn_diode_1d.deck
```

This writes `pn_diode_1d.csv`, `.nitsweep` and `.nitdev` next to where you ran it. To open the GUI:

```bat
build\release\nitcad-app.exe
```

Open a deck from `examples\` and press **F5** to run it.

## 6. Optional: PARDISO solver (Intel MKL)

Large 3D decks (for example `decks\finfet_3d.deck`) should use PARDISO. NiTCAD never links MKL; it loads
`mkl_rt.3.dll` at run time.

1. Download the NuGet packages `intelmkl.redist.win-x64` and `intelopenmp.redist.win` from nuget.org (a `.nupkg`
   is a zip file).
2. Copy every DLL from `runtimes\win-x64\native\` in both packages into one folder, for example `C:\mkl`.
   `mkl_rt.3.dll` and `libiomp5md.dll` must both be there.
3. To run the PARDISO tests, reconfigure with the path:

   ```bat
   cmake --preset release -DNITCAD_MKL_RUNTIME=C:\mkl\mkl_rt.3.dll
   ctest --preset release
   ```

   Without this option the PARDISO tests are skipped.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `VCPKG_ROOT is not set` | Run `set VCPKG_ROOT=C:\vcpkg` in the same prompt before `cmake --preset`. |
| `VCPKG_ROOT ... is not a git checkout at the pinned commit` | The prompt is using Visual Studio's vcpkg, or your checkout is at a different commit. Reset `VCPKG_ROOT` and re-run the `git checkout` from step 2. |
| `MSVC ... is older than the verified reference toolset 19.51` | Install MSVC toolset 14.51 from the Visual Studio Installer, or call `vcvarsall.bat x64 -vcvars_ver=14.51`. |
| `NiTCAD is Windows-only` / `requires MSVC` | You are not in an MSVC x64 environment. Use the x64 Native Tools prompt. |
| `Unknown CMake command` or preset not found | CMake is older than 4.3. Run `cmake --version` and upgrade. |
