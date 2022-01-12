@ECHO OFF
set PATH=C:\OpenSees-FPA\Win64\bin;C:\OpenSees-Solvers\tcl\bin;C:\OpenSees-Solvers\hdf5;C:\OpenSees-Solvers\bin;C:\OpenSees-Solvers\plugins;%PATH%
echo source main.tcl | OpenSees.exe %1
pause