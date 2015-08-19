@echo off
set cygwin64dir=C:\cygwin64

echo Directories:
echo %cygwin64dir%
echo.
echo Make sure that the folders are set ok.

set path=%cygwin64dir%\bin;%path%
set CFLAGS=-lm 
REM -g -O0 -Wall -Wextra
DEPS = ../src/wrapwindfield.c ../src/zephyros_config.c ../src/util.c ../src/fields.c ../src/turbulence.c ../src/interpolation.c ../src/particles.c ../src/coordinates.c ../src/edr.c ../src/func.c ../src/ltqnorm.c ../src/digamma.c

@echo on
gcc -o wrapwindfield_matlab.exe %DEPS% %CFLAGS%
@echo off

echo.
echo.
PAUSE
