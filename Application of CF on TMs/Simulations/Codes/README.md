# Instructions to run "TMS_actual_simulations.R" script

1)  Make sure to use the correct working directory (wd) where all the files from this folder are present (**Change Line 50** in "TMS_actual_simulations.R" script).
2)  Make sure to use the correct working directory (curr_wd) where all the files from this folder are present. Change the lines corresponding to **cur_wd** in "TM_1_simulation.m", "TM_2_simulation.m", "TM_3_simulation.m", "TM_1_odes.m", "TM_2_odes.m", and "TM_3_odes.m" scripts to reflect this change.
   
3)  The aruguments for the "Simulate_TM_1" function in "TMS_actual_simulations.R" are:
   
    a) arg1 - wd
    
    b) arg2 - exchange rates (either 3.2, 320 or 3200 as input)
    
    c) arg3 - WT/Gene KO (either "WT", "B", "E", "X", "A" or "Z" as input)
    
4)  The aruguments for the "Simulate_TM_2" function in "TMS_actual_simulations.R" are:
   
    a) arg1 - wd
    
    b) arg2 - exchange rates (either 3.2, 320 or 3200 as input)
    
    c) arg3 - WT/Gene KO (either "WT", "X", "A" or "Z" as input)

5)  The aruguments for the "Simulate_TM_3" function in "TMS_actual_simulations.R" are:
   
    a) arg1 - wd
    
    b) arg2 - exchange rates (either 3.2, 320 or 3200 as input)
    
    c) arg3 - WT/Gene KO (either "WT", "X" or "A" as input)

6) Run all these cases one at a time to generate the desired output files.


