@echo off
set cygwin64dir=C:\cygwin64

echo Directories:
echo %cygwin64dir%
echo.
echo Make sure that the folders are set ok.

set path=%cygwin64dir%\bin;%path%

@echo on
wrapwindfield_matlab.exe
@echo off

echo.
echo.
