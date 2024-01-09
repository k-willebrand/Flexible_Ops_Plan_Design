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
   2. Generates monthly precipitation and temperature time series for each climate state using a K-nn approach (Rajagopalan & Lall, 1999; Skerker et al., 2023)
   3. Derives monthly inflow time series of inflow into the reservoir using a lumped conceptual rainfall-runoff model (Strzepek & Mccluskey, 2010)
   2. **Key folders:** `BMA_code`, `CLIRUN`,`Streamflow_yield_analysis`, `SDP_expansion`  
   
   
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

This model uses the Mwache Dam in Kenya as a case study.


## Key Model and Data Descriptions

### **Dynamic climate uncertainty model**

1. `multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m`: This is the primary script that defines the infrastructure planning
   SDP. This script includes steps that calculate the climate change transition matrix using Bayesian Model Averaging (BMA)
   (see Fletcher et al., 2019) and temperature, precipitation, and runoff monthly time series for each climate state 
   (Rajagopalan & Lall, 1999; Skerker et al., 2023; Strzepek & Mccluskey, 2010).

2. `BMA_results_RCP85_2020-11-14.mat`: Matlab data file containing the BMA results that are used to develop the transition probability matrix for precipitation (see Fletcher et al., 2019).


3. `Mombasa_TandP.mat`: Matlab data file containing monthly time series of temperature and precipitation from 21 modeled General Circulation Models (GCMs) for 1900-2100.


4. `T_Temp_Precip_RCP85B.mat`: Matlab data file containing the calculated precipitation state transition probability matrix. Temperature transition probabilities are assumed to be deterministic. This file is referenced in the infrastructure planning SDP.


5. `runoff_by_state_02Nov2021.mat`: Matlab data file containing the resulting synthetic 200-year-long monthly time series of temperature, precipitation, and runoff across 21 modeled GCMs.


6. `29Sept2021_CalibratorParams.mat`: Calibration parameters for running the CLIRUN model (see Skerker et al., 2023).


### Infastructure planning SDP

1. `multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m`:
   
   This is the primary script that defines the infrastructure planning SDP. 
   The parameter setup for the file defines parameter and inputs defined in `Multiflex_Wrapper_Nov2021_cprime.m`. 
   The run parameters specify whether to run different components of the model. 


2. `Multiflex_Wrapper_Nov2021_cprime.m`:

   This is the primary script required for interacting with the infrastructure planning SDP. 
   The script is a wrapper that iteratively runs `multiflex_sdp_climate_StaticFlex_DetT_Nov2021`. 
   For one set of cost model parameters (c' and discount rate can be an array), the script does the following:

      1. Finds the optimal flexible design dam size (with static and flexible operations) 
      2. Finds the optimal static dam size (with static and flexible operations) 
      3. Finds the optimal flexible planned dam (with static and flexible operations) 
      4. Runs and saves the forward simulation results for each of the with the six flexible/static infrastructure alternatives


### Reservoir operation SDP

1. `cluster_parfor_main_script_SDP_Evap_SteadyState.m`

   This is the main script utilized to optimize and simulate adaptive (flexible) and non-adaptive (static) reservoir 
   operations using SDP based on Giuliani et al. (2016). This script is parallized. Multiple jobs can be submitted to high performance computing clusters 
   (e.g., across different dam capacities to effeciently optimize and simulate operating policies across 
   different discrete climate states). Postprocessing is then performed on intermediate simulated reservoir operating 
   files to save results in a file format compatible with the infrastructure planning SDP.


2. `cluster_main_script_SDP_postProcessing.m`

   A postprocessing script that reformats and saves the expected reservoir operations shortage costs from
   `cluster_parfor_main_script_SDP_Evap_SteadyState.m` into an file naming convention and structure for use with the 
   infrastructure planning SDP.


3. `Results/Results_SDP_reservoir_ops`
   
   This folder contains previously calculated expected shortage costs for 
   each candidate dam capacity under adaptive (flexible) and non-adaptive (static) operations. These files are used when
   running the infrastructure planning SDP with previously calcuated results from the reservoir operation SDP (recommended).

## Steps to Run the Model

For specific details on how to interact with the model, please refer to the detailed script and parameter descriptions 
provided at the start of key scripts (e.g., `multiflex_sdp_climate_StaticFlex_DetT_Nov2021`) 

To run this model to recreate the results from Willebrand et al. (2024):

1. Download the Github repository


2. Add the folder to your local path. This can be acheived in Matlab by running the following line:
```
addpath(genpath('local directory/Flexible_Ops_Plan_Design'))
```


3. Update the default file path names to correspond with your local path


4. Select method for obtaining the climate state transition probability matrix for the infrastructure planning SDP by either:

   1. Running the script `multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m` where the parameters `runRunoff`, `runTPts`, 
   and `runoffPostProcess` should be set to `true` (default: False). This will recalculate and save 
   use with the infrastructure planning SDP.

   2. Using the previously calculated probability matrix Matlab data file where the parameters `runRunoff`, `runTPts`, 
   and `runoffPostProcess` should be set to `false` in `multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m`(recommended to reduce computational time). 


5. Obtain the expected shortage cost results for each candidate dam capacity and infrastructure operating alternative by either:

   1. Re-running the reservoir operations SDP via `cluster_parfor_main_script_SDP_Evap_SteadyState.m`. 
      1. This script is parallized. Multiple jobs can be submitted to high performance computing clusters 
   (e.g., across different dam capacities to effeciently optimize and simulate operating policies across 
   different discrete climate states). 
      2. Postprocess the resulting saved files using `cluster_main_script_SDP_postProcessing.m` for use with the infrastructure planning SDP.
   2. Using the previously saved shortage cost results obtained from running the reservoir operation SDP. 
      1. The infrastructure planning SDP is configured to use preloaded shortage costs data output from the 
         reservoir operations SDP to reduce the computational time required to repeatedly rerun the reservoir operations SDP. 
      2. The use of preloaded data is specified by `runParam.calcShortage = false` in `multiflex_sdp_climate_StaticFlex_DetT_Nov2021.m`. 
      3. Preloaded data is contained in the folder `Results/Results_SDP_reservoir_ops`.


6. Run the infrastructure planning SDP

   1. Update the cost model parameters (c' water shortage cost parameter and discount rate) in `Multiflex_Wrapper_Nov2021_cprime.m`. 
      2. Cost parameters can be specified as arrays to test different cost scenarios. 
      3. This is helpful when conducting senstivity analysis for c' and discount rate. For example, c' sensitivity analysis results 
      from Willebrand et al. (2024) are saved in 'Results/Results_cprime_sensitivity'
   2. To consider alternative dam capacity assumptions (e.g., specify initial dam capacity as 60 MCM), update conditional constraints on dam sizes in `Multiflex_Wrapper_Nov2021_cprime.m`. 
   3. Run `Multiflex_Wrapper_Nov2021_cprime.m`. Results by default will be saved in `Results/Results_SDP_expansion`.

7. Analyze the results

   1. Final results are saved in the `Results` folder for analysis
   2. Consider using the included plotting utilities in the `Plots` folder


## References
Fletcher, S., Lickley, M., & Strzepek, K. (2019). Learning about climate change uncertainty enables flexible water infrastructure planning. *Nature Communications*, 10(1), Article 1. https://doi.org/10.1038/s41467-019-09677-x


Giuliani, M., Y. Li, A. Cominola, S. Denaro, E. Mason, and A. Castelletti (2016), A Matlab toolbox for designing Multi-Objective Optimal Operations of water reservoir systems, *Environmental Modelling & Software*, 85, 293–298


Rajagopalan, B., & Lall, U. (1999). A k-nearest-neighbor simulator for daily precipitation and other weather variables. *Water Resources Research*, 35(10), 3089–3101. https://doi.org/10.1029/1999WR900028


Skerker, J. B., Zaniolo, M., Willebrand, K., Lickley, M., & Fletcher, S. M. (2023). Quantifying the Value of Learning for Flexible Water Infrastructure Planning. *Water Resources Research*, 59(6), e2022WR034412. https://doi.org/10.1029/2022WR034412


Strzepek, K. M., & Mccluskey, A. L. (2010). Modeling the Impact of Climate Change on Global Hydrology and Water Availability. World Bank Group. https://assets.publishing.service.gov.uk/media/57a08b18ed915d3cfd000b20/EACC-Hydrology.pdf


Willebrand, K., Zaniolo, M., Skerker, J.B., and Fletcher, S.M., Assessing the interacting value of flexible planning, design, and operations in water supply infrastructure, *Water Resources Research*, 2024. [In Review]


