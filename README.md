# SAS-N <a href="https://github.com/tamnva/SAS_N/sas_logo.svg"><img src="sas_logo.svg" align="right" height="120" /></a>

[![DOI](https://zenodo.org/badge/509562497.svg)](https://zenodo.org/badge/latestdoi/509562497)
- SAS-N is a catchment scale nitrogen transport model using the storage age selection function.
- This repository contains the source code (./src) and test case (./test) from the SAS_N model. 
- The output folder (e.g., ./test/output) contains the following files.

  ***Nbalance.txt***

     - Columns are years (from 1950-2014)

     - Rows 1 to 8 are N surplus (kg/ha/year), soil N store  (kg/ha), total soil N removal (soil N denitrification + N leaching to the groundwater)  (kg/ha), soil N denitrification (kg/ha/year), soil N leaching (kg/ha/year), riverine N export (kg/ha/year), groundwater N storage (kg/ha), and groundwater denitrification (kg/ha/year), respectively. This repeats 30 times (from the best 30 simulations), in other words, the first 8 rows contain results from the first simulation, the next 8 rows contain results from the next simulation, and so on.

       
  
  ***RMSE*** 
  
     - The model performance, the root means square error (mg/L), of the 30 best simulations (row 1 to 30).
  
       
  
  ***Sim_C*** 
  
     - Columns are the years from 1950-2014
  
     - Rows 1 to 30 contain the simulated nitrate (N-NO3) concentrations from the subsurface (mg/L) from the 30 best simulations.
  
       
  
  ***TTD.txt*** 
  
     - The mean transit times (year) from the best 30 simulations (rows 1 to 30).
  


​     
