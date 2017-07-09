@echo off
%~dp0luajit %~dp0lua\trepl\th -e "os.execute('start http://localhost:8000');require('display.start')"
exit /b %ERRORLEVEL%
