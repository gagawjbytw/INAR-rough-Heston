@echo off
echo Creating build directory if it doesn't exist...
if not exist "..\bin" mkdir "..\bin"

echo Compiling monte_carlo_inar_multi_strikes.cpp...
g++ ..\src\monte_carlo_inar_multi_strikes.cpp -o ..\bin\monte_carlo_inar_multi_strikes.exe -std=c++11 -pthread -O2
if %errorlevel% neq 0 (
    echo Compilation failed!
    pause
    exit /b 1
)
echo Compilation successful!
echo.
echo Running monte_carlo_inar_multi_strikes.exe...
echo.
cd ..\bin
monte_carlo_inar_multi_strikes.exe
echo.
pause 