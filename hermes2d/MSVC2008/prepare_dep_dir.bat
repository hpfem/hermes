	echo off
	echo Creating a directory structure for dependecy libraries
	set PDD_OLD_DIR=%CD%

	if exist inttypes.h goto INTTP_FOUND
	echo A file 'inttypes.h' not found. Script not executed in MSVC2008 directory. Please, copy contents of a directory 'MSVC2008' to hermes2d and run the script 'prepare_dep_dir.bat'.
	goto END
:INTTP_FOUND

	echo This script has to be executed in MSVC2008 directory.
	cd ..\..
	if not exist dependecies goto CREATE_DIR
	echo Directory structure seems to exist already. Quiting.
	goto END
	
:CREATE_DIR
	mkdir dependecies
	mkdir dependecies\bin
	mkdir dependecies\include
	mkdir dependecies\lib
	cd dependecies\include
	copy %PDD_OLD_DIR%\inttypes.h %CD%
	cd ..\bin
	echo Please add following to you PATH variable: %CD%

:END
	cd %PDD_OLD_DIR%


