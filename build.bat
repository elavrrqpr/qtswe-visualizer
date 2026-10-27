@echo off
REM Windows build (needs gcc from MinGW-w64 / MSYS2 on PATH)
gcc -O3 -std=gnu11 -fopenmp -D__USE_MINGW_ANSI_STDIO=1 -o qtswe.exe main.c qtswe.c -lws2_32 -lm
if errorlevel 1 (echo BUILD FAILED) else (echo OK: qtswe.exe built)
