@echo off
cd ..
echo Compiling monte_carlo_inar_all_options.cpp...
g++ -std=c++17 src/monte_carlo_inar_all_options.cpp -o monte_carlo_inar_all_options.exe

if %ERRORLEVEL% NEQ 0 (
    echo Compilation failed!
    echo Please make sure you have g++ installed and added to your PATH
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Running the program...
echo.
monte_carlo_inar_all_options.exe

echo.
pause 
