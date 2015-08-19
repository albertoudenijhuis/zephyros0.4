@echo off
set cygwin64dir=C:\cygwin64

echo Directories:
echo %cygwin64dir%
echo.
echo Make sure that the folders are set ok.

set path=%cygwin64dir%\bin;%path%
set CFLAGS=-lm 
REM -g -O0 -Wall -Wextra
set DEPS=../src/wrapradarfilter_matlab.c ../src/wrapradarfilter.c ../src/radarfilter.c ../src/fields.c ../src/turbulence.c ../src/edr.c ../src/func.c ../src/interpolation.c ../src/coordinates.c ../src/digamma.c ../src/ltqnorm.c ../src/util.c ../src/zephyros_config.c

@echo on
gcc -o wrapradarfilter_matlab.exe %DEPS% %CFLAGS%
@echo off

echo.
echo.
PAUSE
