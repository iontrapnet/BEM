@echo off
%~dp0luajit -e "os.execute('start http://localhost:8000');arg={...};require('display.start')"
exit /b %ERRORLEVEL%
