@echo off
rem Builds build\installer\NiTCAD-<version>-setup.exe: Release nitcad and nitcad-app, staged with
rem `cmake --install`, packed by Inno Setup 6 (winget install JRSoftware.InnoSetup).
rem Run from a "x64 Native Tools" prompt with VCPKG_ROOT set, or after vcvars64.bat.
setlocal
cd /d "%~dp0.."
cmake --preset release || exit /b 1
cmake --build --preset release -j 6 --target nitcad nitcad-app || exit /b 2
if exist build\stage rmdir /s /q build\stage
cmake --install build\release --prefix build\stage || exit /b 3
for /f "tokens=4 delims=( " %%v in ('findstr /b /c:"project(NiTCAD VERSION" CMakeLists.txt') do set VERSION=%%v
set ISCC=%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe
if not exist "%ISCC%" set ISCC=%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe
if not exist "%ISCC%" (echo Inno Setup 6 not found: winget install JRSoftware.InnoSetup & exit /b 4)
"%ISCC%" /Q /DAppVersion=%VERSION% packaging\nitcad.iss || exit /b 5
echo Built build\installer\NiTCAD-%VERSION%-setup.exe
