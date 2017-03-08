:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Monteverdi launcher to set up proper environment
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
@echo off
set CURRENT_SCRIPT_DIR=%~dp0
setlocal
call %CURRENT_SCRIPT_DIR%otbenv.cmd

:: Set current dir to HOME dir because Monteverdi generates temporary files and need write access
cd %HOMEDRIVE%%HOMEPATH%

:: Start Monteverdi
start "Monteverdi Viewer" /B "%CURRENT_SCRIPT_DIR%bin\monteverdi.exe" %*
endlocal
