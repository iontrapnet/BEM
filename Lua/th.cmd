@echo off
%~dp0luajit %~dp0lua\trepl\th %*
exit /b %ERRORLEVEL%
