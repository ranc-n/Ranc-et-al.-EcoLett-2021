# Ranc-et-al.-EcoLett-2021


Publicly available dataset and program code for:

**N. Ranc, F. Cagnacci & P.R. Moorcroft. 2021. Memory drives the formation of animal home ranges: evidence from a reintroduction. *Ecology Letters* Vol. xx xx-xx. DOI: XXXXXXX
**
<br> 
<br>
The mechanistic movement model is written as a C++ program. The code made available here (1) calculates the likelihood of a specific parameter set given the empirical movement data or (2) simulates movements according to the specific parameter set.\
<br>
<br>

### 1. Inputs
1. *Aspromonte_roedeer_traj.txt* (! published separately!): input files for the model fitting (movement trajectories)
1. *Simulation_traj files*: input files for the movement simulations (release locations for each individual)
1. *Environmental_rasters files*: environmental rasters used for model fitting or movement simulations
1. *Configuration files*: input parameter values

#### Aspromonte_roedeer_traj.txt
The file can be accessed at: 10.5281/zenodo.5189835. This file stores the 6-hour movement trajectory for each animal and is used for model fitting. It includes the following columns:
* *animals_id*: unique animal identifier
* *x*: longitude (epsg 3035; relative to the origin of the environmental rasters; missing data indicated by -9999)
* *y*: latitude (epsg 3035; relative to the origin of the environmental rasters; missing data indicated by -9999)

#### Simulation_traj files
The *Simulation_traj* folder contains two files storing the released locations of each animal. These files are used for the movement simulations and include the same fields as *Ranc_et_al_2021_Aspromonte_roedeer_traj.txt*.
* *Aspromonte_roedeer_release_1y.txt*: released locations for roe deer followed for one year or less.
* *Aspromonte_roedeer_release_2y.txt*: released locations for roe deer followed for more than a year.

#### Environmental_rasters files
The *Environmental_rasters* folder contains all the raster files (asci files, 25m spatial grain) used for model fitting and movement simulations:
* *slope.asc*: slope (standardized).
* *slope_sq.asc*: slope squared.
* *tcd_325grain.asc*: tree cover density (focal average at 325m-grain) (standardized).
* *tcd_325grain_sq.asc*: tree cover density (focal average at 325m-grain) squared.
* *landcover_agri.asc*: areas dominated by agriculture, pasture and anthropogenic areas.
* *landcover_5322.asc*: areas reforested with decidous trees.

#### Configuration files
The *Configuration_files* folder stores the configuration files required to reproduce the results of the fitted models (the parameter values in the config files have been estimated through maximum likelihood) and to simulate movements according to the estimated parameter values:
* *config_best_Mmem_fitting.txt*: fitting of the memory-based model
* *config_best_Mres_fitting.txt*: fitting of the resource-only model
* *config_best_Mmem_simulation_1y.txt*: movement simulations of the memory-based model (individuals monitored for a year or less) 
* *config_best_Mmem_simulation_2y.txt*: movement simulations of the memory-based model (individuals monitored for more than a year) 
* *config_best_Mres_simulation_1y.txt*: movement simulations of the resource-only model (individuals monitored for a year or less) 
* *config_best_Mres_simulation_2y.txt*: movement simulations of the resource-only model (individuals monitored for more than a year) 

The configuration files list 19 parameters:
* *output_directory*: /Path/to/outputs/
* *raster_directory*: /Path/to/Environmental_rasters/
* *trajectory_path*: /Path/to/traj.txt
* *writing_outputs*: true (detailed outputs written); false
* *model_type*: kernel (fitting process), simulation (movement simulation)
* *threshold_approx_kernel*: distance in meters beyond which movement probability = 0 (7000)
* *threshold_memory_kernel*: distance in meters beyond which no learning occurs (1000)
* *resource_names*: vector of environmental layer names
* *selection_coef*: vector of environmental layer selection coefficients
* *r_l*: learning rate of reference memory
* *w_l*: learning rate of working memory
* *r_d*: decay rate of reference memory
* *w_d*: decay rate of working memory
* *r_dist*: spatial scale of learning of reference memory
* *w_dist*: spatial scale of learning of reference memory
* *step_length_dist*: Weibull rate parameter for step length
* *step_length_shape*: Weibull shape parameter for step length
* *n_simulated_steps*: number of simulated steps per trajectory (if model_type=simulation; ignored otherwise)
* *n_simulated_runs*: number of simulated trajectories per individual (if model_type=simulation; ignored otherwise)
<br>


### 2. Outputs
In fitting mode, the program has six outputs:
* *objective_function.csv*: log-likelihood of the model parameter set given the data (single output when writing_outputs=false)
* *objective_function_detail.csv*: step-level log-likelihood of the model parameter set given the data (*r_patch* and *c_patch* translate continuous coordinates into row and column raster coordinates).
* *arenaLogFile.txt*: log file for the creation of the spatial landscape (arena)
* *mainLogFile.txt*: log file for the fitting of the model to all the individuals
* *global_resource.asc*: combined resource selection function based on all the environmental resource layers and associated selection coefficients
* *trajectory_covariates.csv*: intersection of resource selection function and movement trajectory
<br>

In movement simulation mode, the program has xx outputs:
* *simulations.csv*: movement simulations (*run* = 0 for the release locations, 1:n for the simulation run of each individual)
* *mainLogFile.txt*: log file for the fitting of the model to all the individuals
<br>


### 3. Running instructions
To run the code, the program must first be compiled. For example, using the gcc compiler with an unix shell command line:

*cd /Path/to/Code*
*g++ -std=c++11 -O3 Configuration.cpp Dist_lookup_table.cpp Import_traj.cpp Launch_arena.cpp Likelihood.cpp Main.cpp Patch_dynamics.cpp Writing_outputs.cpp -o redistribution_kernel_roedeer*
<br>

The program can then be executed with a given configuration file; for example:
*./redistribution_kernel_roedeer -config /Path/to/config_best_Mmem_fitting.txt*
<br>

For parameter estimation, the program can be called by any optimization algorithm (e.g., MCMC, Particle Swarm Optimization, etc.).
