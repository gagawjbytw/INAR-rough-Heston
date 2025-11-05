@echo off
setlocal

:: Get the directory where the script is located
set "SCRIPT_DIR=%~dp0"
:: Remove trailing backslash
set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
:: Get the project root directory (parent of scripts directory)
for %%I in ("%SCRIPT_DIR%\..") do set "PROJECT_ROOT=%%~fI"

echo Compiling Heston Model Program...

:: Change to project root directory
cd /d "%PROJECT_ROOT%"

:: Create bin directory if it doesn't exist
if not exist "bin" mkdir bin

:: Compile the C++ program
g++ -std=c++17 -O3 -pthread src/classic_Heston_EM.cpp -o bin/heston_model.exe

:: Check if compilation was successful
if %ERRORLEVEL% EQU 0 (
    echo Compilation successful!
    echo Running the program...
    echo.
    bin\heston_model.exe
) else (
    echo Compilation failed!
    pause
    exit /b 1
)

echo.
echo Program execution completed.
pause 