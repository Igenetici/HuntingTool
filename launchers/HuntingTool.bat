@echo off
REM Windows: double-click to open the HuntingTool window.
cd /d "%~dp0.."
where huntingtool-gui >nul 2>nul && (start "" huntingtool-gui & exit /b)
where pyw >nul 2>nul && (start "" pyw huntingtool.py gui & exit /b)
start "" pythonw huntingtool.py gui
