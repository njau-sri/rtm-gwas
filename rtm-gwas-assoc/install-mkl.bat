@echo off

if not defined RTM_GWAS_VERSION (
    set RTM_GWAS_VERSION=unknown
)

set /p my_target=Specify Target Platform [win32/win64]:

rem Distributable Code for Visual Studio 2019
rem https://docs.microsoft.com/en-us/visualstudio/releases/2019/redistribution#visual-c-runtime-files
set "VisualStudioFolder=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community"

set my_vcvar_bat=

if /i "%my_target%" == "win32" (
    rmdir /s /q Release rtm-gwas-assoc-win32
    set "my_vcvar_bat=%VisualStudioFolder%\VC\Auxiliary\Build\vcvars32.bat"
    set "my_vcomp_dll=%VisualStudioFolder%\VC\Redist\MSVC\14.21.27702\x86\Microsoft.VC142.OPENMP\vcomp140.dll"
)

if /i "%my_target%" == "win64" (
    rmdir /s /q x64 rtm-gwas-assoc-win64
    set "my_vcvar_bat=%VisualStudioFolder%\VC\Auxiliary\Build\vcvars64.bat"
    set "my_vcomp_dll=%VisualStudioFolder%\VC\Redist\MSVC\14.21.27702\x64\Microsoft.VC142.OPENMP\vcomp140.dll"
)

if "%my_vcvar_bat%" == "" (
    echo ERROR: invalid target platform: %my_target%
    pause
    exit /b 1
)

call "%my_vcvar_bat%"

MSBuild.exe -m -p:Configuration=Release

if errorlevel 1 (
    pause
    exit /b 1
)

if "%my_target%" == "win32" (
    mkdir rtm-gwas-assoc-win32
    copy "%my_vcomp_dll%" rtm-gwas-assoc-win32\
    copy Release\rtm-gwas-assoc.exe rtm-gwas-assoc-win32\
)

if "%my_target%" == "win64" (
    mkdir rtm-gwas-assoc-win64
    copy "%my_vcomp_dll%" rtm-gwas-assoc-win64\
    copy x64\Release\rtm-gwas-assoc.exe rtm-gwas-assoc-win64\
)

pause

set my_target=
set my_vcvar_bat=
set my_vcomp_dll=
set VisualStudioFolder=
