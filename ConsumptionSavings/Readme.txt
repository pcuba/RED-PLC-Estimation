# README
# This file details instructions to reproduce the solution of the
# Consumption-Savings Model with Borrowing Constraint from Appendix F.
#
# "Piecewise-Linear Approximations and Filtering
# for DSGE Models with Occasionally Binding Constraints," by Aruoba,
# Cuba-Borda, Higa-Flores, Schorfheide and Villavazo (2020).
#
#===============================================================================

A. Folder structure
-------------------
The scripts associated for each solution method are in the subdirectories  FiPIt,
OccBin, PLC. Auxiliary codes are in the subdirectory External.

B. OccBin Installation
----------------------
You need to download OccBin from https://www.matteoiacoviello.com/research.htm.
Note that OccBin might crash for versions of Dynaer >4.5. See OccBin/setpathdynare4.m

C. MAIN SCRIPT
--------------
The file script_solve_borrcon.m solves the model using FiPIt, OccBin and PLC. It
also produces decision rules for each solution and simulated paths.
