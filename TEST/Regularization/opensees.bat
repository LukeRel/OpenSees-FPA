@ECHO OFF
set PATH=C:\OpenSees\Win64\bin;C:\OpenSees-Solvers\tcl\bin;C:\OpenSees-Solvers\hdf5;C:\OpenSees-Solvers\bin;C:\OpenSees-Solvers\plugins;%PATH%
OpenSees.exe %1
pause