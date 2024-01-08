# Flexible_Ops_Plan_Design Repository 

This repository presents a Matlab implementation that models three types of flexibility in water supply infrastructure 
projects together: planning, design, and operations. The approach integrates a dynamic model of climate uncertainty with 
joint optimization of stochastic dynamic control policies for infrastructure planning and operation. This allows for
an assessment of the interactions and the relative role different types of flexibility can play in mitigating cost and 
water supply reliability risk under climate change uncertainty.

The code in this repository replicates the results in:

* Willebrand, K., Zaniolo, M., Skerker, J.B., and Fletcher, S.M., Assessing the interacting value of flexible planning, design, and operations in water supply infrastructure, Water Resources Research, 2024.

Code documentation by: Keani Willebrand, Stanford University. Last updated: January 2024

## Model Overview

Central to this repository is the development of a nested stochastic dynamic programming (SDP) model that jointly optimizes
dam infrastructure planning and operating policies under dynamic climate uncertainty. 

This repository is comprised of three primary modeling components and a plotting utility for analysis:

### 1. **Dynamic climate uncertainty model**
   1. Calculates and saves the precipitation state transition probability matrix that is referenced in the infrastructure planning SDP
   2. Generates monthly precipitation and temperature time series for each climate staet using a K-nn approach 
   3. Derives monthly inflow time series of inflow into the reservoir using the a lumped conceptual rainfall-runoff model
   2. **Key folders:** `BMA_code`, `CLIRUN`,`Streamflow_yield_analysis`  
   
   
### 2. **Reservoir operation SDP**
   1. Derives optimal adaptive (flexible) and non-adaptive (static) operating policies for each candidate dam design and climate state
   2. Simulates optimal operating policies to calculate expected shortage costs for each reservoir capacity and climate state
   3. Saves the expected shortage costs in files that are referenced in the infrastructure planning SDP
   2. **Key folder:** `SDP_reservoir_ops`
   

### 3. **Infrastructure planning SDP**
   1. Central component of the model
   2. Derives optimal policies that determine initial infrastructure and subsequent capacity expansion decisions as a function of the current state of the system.
   3. Simulates optimal infrastructure planning policies to calculate expected total costs by end of century 
   4. Saves the expected shortage costs and simulated actions in files for analysis and plotting 
   5. **Key folder:** `SDP_expansion`


### 4. **Plotting and analysis**  
   1. Plotting utilties for recreating the figures published in Willebrand et al. (2024)
   2. **Key folder:** `Plots`

This model uses the Mwache Dam in Kenya as a case study example.


## Key Model and Data Descriptions

### **Dynamic climate uncertainty model** 


2. `T_Temp_Precip_RCP85B.mat`: Matlab data file containing the calculated precipiation state transition probability matrix. This file is referenced in the infrastructure planning SDP.


3. `runoff_by_state_02Nov2021.mat`: Matlab data file containing the resulting 200-year-long monthly time series of temperature, precipitation, and runoff across 21 modeled General Circulation Models (GCMs).

### Infastructure planning SDP

1. `multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m`:
   
   This is the primary script that defines the infrastructure planning SDP. 
   The parameter setup for the file explicitly defines some parameters and inputs others from Multiflex_OptimizeSimulate.m. The run parameters at the top of the script are used to specify which parts of the script to run different components of the model. 


2. `Multiflex_Wrapper_Nov2021_cprime.m`:

   This is the primary script required for interacting with the infrastructure planning SDP. 
   The script is a wrapper that iteratively runs `multiflex_sdp_climate_StaticFlex_DetT_Nov2021`. 
   For one set of cost model parameters (c' and discount rate can be an array), the script does the following:

      1. Find the optimal flexible design dam size (with static and flexible operations) 
      2. Find the optimal static dam size (with static and flexible operations) 
      3. Find the optimal flexible planned dam (with static and flexible operations) 
      4. Run and save the forward simulation results for each of the with the six flexible/static infrastructure alternatives



### Reservoir operation SDP

1. `cluster_parfor_main_script_SDP_Evap_SteadyState.m`

   This is the main script utilized to optimize and simulate adaptive (flexible) and non-adaptive (static) reservoir 
   operations using SDP. This script is parallized. Multiple jobs can be submitted to high performance computing clusters 
   (e.g., across different dam capacities to effeciently optimize and simulate operating policies across 
   different discrete climate states). 


   Then, postprocessing is performed to save results in a format compatible with the infrastructure planning SDP code.

   

2. `cluster_main_script_SDP_postProcessing.m`

   A postprocessing script that reformats and saves the expected reservoir operations shortage costs from
   `cluster_parfor_main_script_SDP_Evap_SteadyState.m` into an appropriate matrix structure for use with the 
   infrastructure planning SDP.

## Steps to Run the Model

For specific details on how to interact with the model, please refer to the detailed script and parameter descriptions 
provided at the start of key scripts (e.g., `multiflex_sdp_climate_StaticFlex_DetT_Nov2021`) 

To run this model to recreate the results from Willebrand et al. (2024):

1. Download the Github repository
2. Add the folder to your local path. This can be acheived in Matlab by running the following line:
```
addpath(genpath('local directory/Flexible_Ops_Plan_Design'))
```
3. Update the default file path names to correspond with your updated local path
4. Obtain the climate change transition probability matrix for the infrastructure planning SDP by either:
   1. 
5. Obtain the expected shortage cost results for each candidate dam capacity and infrastructure operating alternative (recommended use of parallelization)
   1. use the previously saved shortage cost results obtained from the simulatation of optimized adaptive (flexible) and non-adaptive (static) reservoir operation policies.
   2. Run 
6. Run the infrastructure planning SDP
   1. Select
