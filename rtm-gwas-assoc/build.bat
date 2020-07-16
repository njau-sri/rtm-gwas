@echo off

set /p MyTarget=Specify Target Platform [win32/win64]:

rem Distributable Code for Visual Studio 2019
rem https://docs.microsoft.com/en-us/visualstudio/releases/2019/redistribution#visual-c-runtime-files
set "VisualStudioFolder=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community"
set "RedistVersion=14.26.28720"

set MyVcVarBat=

if /i "%MyTarget%" == "win32" (
    rmdir /s /q Release
    set "MyVcVarBat=%VisualStudioFolder%\VC\Auxiliary\Build\vcvars32.bat"
    set "MyVcOmpDll=%VisualStudioFolder%\VC\Redist\MSVC\%RedistVersion%\x86\Microsoft.VC142.OPENMP\vcomp140.dll"
)

if /i "%MyTarget%" == "win64" (
    rmdir /s /q x64
    set "MyVcVarBat=%VisualStudioFolder%\VC\Auxiliary\Build\vcvars64.bat"
    set "MyVcOmpDll=%VisualStudioFolder%\VC\Redist\MSVC\%RedistVersion%\x64\Microsoft.VC142.OPENMP\vcomp140.dll"
)

if "%MyVcVarBat%" == "" (
    echo ERROR: invalid target platform: %MyTarget%
    pause
    exit /b 1
)

call "%MyVcVarBat%"

for /f "delims=" %%x in (../VERSION) do set RTM_GWAS_VERSION=%%x

MSBuild.exe -m -p:Configuration=Release

if errorlevel 1 (
    pause
    exit /b 1
)

if "%MyTarget%" == "win32" (
    mkdir Release\bin
    copy "%MyVcOmpDll%" Release\bin\
    copy Release\rtm-gwas-assoc.exe Release\bin\
)

if "%MyTarget%" == "win64" (
    mkdir x64\Release\bin
    copy "%MyVcOmpDll%" x64\Release\bin
    copy x64\Release\rtm-gwas-assoc.exe x64\Release\bin
)

pause

set MyTarget=
set MyVcVarBat=
set MyVcOmpDll=
set RedistVersion=
set VisualStudioFolder=
