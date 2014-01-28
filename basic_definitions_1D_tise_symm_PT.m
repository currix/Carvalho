##   basic_definitions.m
##
##   Basic definitions and launcher for the set of programs solving the 1D TISE for a symmetric potential
##     
## by Currix TM
##
## Uncomment to Clear figure
## clf
##
##
## Verbosity flag
global iprint = 1;
##
##
## Physical constants
global hbarc = 197.32858; # MeV fm
global amu = 938.92635;   # MeV / c^2
global hsqoamu = 41.4713768; # MeV fm^2 
##
##
## Define global variables characterizing the 1D system
## Spatial grid
global xmin = -35; # (fm)
global xmax = 35;  # (fm)
global npoints = 1003; 
global xgrid  = linspace(xmin,xmax,npoints); # Interval comprising ends with npoints points (fm)
global x_step = (xmax-xmin)/(npoints-1); # (fm)
##
## Reduced mass
global red_mass = 0.947 # Reduced mass in amu
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf("xmin = %f fm,\t xmax = %f fm,\t n_points = %d,\t h = %f fm\n", xmin, xmax, npoints, x_step)
  printf(" reduced mass = %f amu\n", red_mass);
endif
## 
##
## Momentum grid (fm-1)
global k_min = 0.025; # fm-1
global k_max = 2.5; # fm-1
global n_k_points = 100;
global k_values = linspace(k_min, k_max, n_k_points); ## Vector with k_values (fm-1)
global E_values = (k_values.*hbarc).^2/(2*red_mass*amu) ## Energy values (MeV)
##
## Energy grid (MeV) (uncomment next lines and comment the  k grid lines)
##global E_min = 0.1;  # MeV
##global E_max = 25.0; # MeV
##global n_E_points = 150;
##global E_values = linspace(E_min, E_max, n_E_points); ## Vector with energy values
##global k_values = sqrt((2*red_mass*amu).*E_values)/hbarc; ## Vector with energy values
##
## Define global variables characterizing the 1D potential
##
## Woods-Saxon Potential parameters
global V_pt = 50.0; # Potential Depth (MeV)
global a_pt = 2.0; # Potential Diffusivity (fm)
##
## Potential values
global vpot =  poeschl_teller_1D(xgrid);
##
## Match Point (jeje, bound states calculation)
global match_p = 400;
##
## Output
if ( iprint >= 1 )
  disp("####################################################################");
  printf(" 1D Poeschl-Teller Potential \n");
  printf(" V_pt = %f MeV, a_pt = %f fm \n", V_pt, a_pt);
  printf(" match point = %d, x_M = %f fm\n", match_p, xgrid(match_p))
endif
##
##
## Matrix diagonalization states
global eigenvectors_file = "ho_eigenvectors_N120.dat";
global dim_N = 120;
global bound_states = 1;
global pseudo_states = dim_N - bound_states;
##
##
#################################################################
#################################################### Calculations
#################################################################
##
## Numerov algorith to compute Bound eigenvalues and eigenstates
##
## Save bound states wave function
global iwf_bound_save = 1;
##
##  Bound Eigenstates filenames wf_filename_1.dat wf_filename_2.dat ... 
global wf_filename = "wf_octave_bound";
## 
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("############  Bound states eigensystem  ############################");
  disp("####################################################################");
  disp("");
endif
##
bound_states_eigensystem_Numerov_symm_pot_1D_tise;
##
if ( iprint >= 1 )
  disp(" ");
endif
##
##
## Continuum Eigenstates
##
## Save continuum states wave function
global iwf_cont_save = 1;
##  Continuum Eigenstates filenames 
wf_filename = "wf_octave_continuum";
##
## Save S matrix
global iSM_save = 1; 
##  S Matrix Filename
global smat_filename = "smatrix_octave_symm.dat";
##
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("###################  Continuum states ##############################");
  disp("####################################################################");
  disp("");
endif
##
continuum_symm_states_Numerov_symm_pot_1D_tise;
##
##
## Test sum rules
##
## Total Strength
##
global iSum_Rules_save = 1;
##
## Continuum symmetrized states
wf_filename = "wf_octave_continuum";
##
if ( iprint >= 1 )
  disp("");
  disp("####################################################################");
  disp("############### Bound States Sum Rules #############################");
  disp("####################################################################");
  disp("");
endif
##
bound_states_sum_rules_Numerov_symm_pot_1D_tise;
##
##
## Response function dB/dE computed with continuum states
##
global i_E = 1; ## 1 -> E1  :: 2 -> E2
##
global isave_dBdE = 1;
##
## Continuum symmetrized states
wf_filename = "wf_octave_continuum";
##
## Response function filename
global dBdE_filename = "response_function_continuum";
##
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("############### Response Function (continuum) ######################");
  disp("####################################################################");
  disp(" ");
endif
##
dBdE_pure_cont_symm_states_Numerov_symm_pot_1D_tise;
##
##
## Pseudodensity calculation
##
##
## Save density function
global idensity_save = 1;
##
##  Continuum Eigenstates filenames 
wf_filename = "wf_octave_continuum";
##
##  Quasidensity filename 
global qdensity_filename = "wfc_rho";
##
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("##################### Quasidensity matrix ##########################");
  disp("####################################################################");
  disp(" ");
endif
##
pseudodensity_symm_states_Numerov_symm_pot_1D_tise;
##
##
##
##  
## Response function dB/dE computed with continuum states
##
global isave_pseudo_dBde = 1;
## Save dBdE filenames
global dBdE_pseudo_filename = "dBde_E1_rho_ISQW.dat";
##
i_E = 1; ## 1 -> E1  :: 2 -> E2
##
##  Quasidensity filename 
qdensity_filename = "wfc_rho";
##
## Read Pseudostates transition moment from Fortran code
global response_function_pseudostate_file = "ho_E2_TM_N120_1.dat";
##
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("################# Pseudostate Response Function ####################");
  disp("####################################################################");
  disp(" ");
endif
##
dBdE_pseudostates_symm_states_Numerov_symm_pot_1D_tise
##
##
##  
## Bins computed by the average method
##
global ibins_save = 1;
global wf_bins = "wfc_octave_bin";
##
## Bin definition
##
global bin_energies = [0.227 0.336 0.908 1.217]; ## Energies for which bins would be computed. (MeV)
##
global delta_k = 0.09; ## bin width in fm-1
##
global n_funs_bin = 80; ## number of functions per bin
##
if ( iprint >= 1 )
  disp(" ");
  disp("####################################################################");
  disp("######################## Average Method bins #######################");
  disp("####################################################################");
  disp(" ");
endif
##
bins_continuum_normalized_symmetry_Numerov_symm_pot_1D_tise
