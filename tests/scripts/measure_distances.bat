@echo off
REM Path to the VMD executable
set VMD_PATH="C:\Program Files (x86)\University of Illinois\VMD\vmd.exe"

REM Path to your TCL script
set TCL_SCRIPT="C:\Users\Thomas\Documents\workspace\research\OpenABC_RNA\tests\scripts\measure_distances.tcl"

REM Path to your PDB file
set PDB_FILE="C:\Users\Thomas\Documents\workspace\research\OpenABC_RNA\tests\rna-pdb-files\1hr2_simulation.pdb"

REM Path to your DCD file
set DCD_FILE="C:\Users\Thomas\Documents\workspace\research\OpenABC_RNA\tests\output.dcd"

REM Residues to measure distances between (replace 10 and 50 with your residue IDs)
set RESID1=3
set RESID2=155

REM Run VMD with the TCL script and input files
%VMD_PATH% -e %TCL_SCRIPT% -args %PDB_FILE% %DCD_FILE% %RESID1% %RESID2%

REM Pause to keep the command prompt open after execution
pause
