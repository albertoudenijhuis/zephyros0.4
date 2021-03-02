@echo off
set cygwin64dir=C:\cygwin64

echo Directories:
echo %cygwin64dir%
echo.
echo Make sure that the folders are set ok.

set path=%cygwin64dir%\bin;%path%
set CFLAGS=-lm -lcxsparse -lnlopt-0
REM -g -O0 -Wall -Wextra
set DEPS=../src/wrapwindvectors_matlab.c ../src/wrapwindvectors.c ../src/windvectors_lwm.c ../src/windvectors_fdvar.c ../src/util.c ../src/radarfilter.c ../src/fields.c ../src/edr.c ../src/func.c ../src/interpolation.c ../src/coordinates.c ../src/digamma.c ../src/ltqnorm.c ../src/zephyros_config.c

@echo on
gcc -o wrapwindvectors_matlab.exe %DEPS% %CFLAGS% -L..\winlib -I..\winlib -I%cygwin64dir%\usr\include\suitesparse
@echo off

echo.
echo.
PAUSE
