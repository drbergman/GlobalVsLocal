function p = basePars_FGFR3()

plot_properties.plotFigs = false;               
plot_properties.makeMovie = false;   
plot_properties.plotLocations = true;
plot_properties.plot_every = 1;
plot_properties.plot_offset = 0;

main.dt = (1/6) / 24; % number of days per step
main.cell_width = 20e-6; % in meters; bladder cancer cell is about 20micrometers in diameter
main.track_all_phiD = false; % whether or not to track all phiD values
main.plot_every_pde_step = false;
main.plot_inhibitor = false;

%% neighbor parameters
main.occmax = 20; % below this threshold, a tumor/immune cell can divide; at and above, too many neighbors and so doesn't proliferate
main.neighbors = allCombos(-1:1,-1:1,-1:1,'matlab');
main.neighbors(all(main.neighbors==0,2),:) = []; % don't count self as neighbor

%% tumor parameters
main.alpha1 = 1; % proliferation rate of tumor cells

main.alpha2 = .4; % (twice the) max increase to prolif probability for tumor cells with FGFR3 signaling; twice because 0<=phiD=D_A/RT<=0.5

main.min_prolif_wait = 9/24; % number of days all cells must wait at minimum between proliferations; comes from Patient-calibrated agent-based modelling of ductal carcinoma in situ (DCIS): From microscopic measurements to macroscopic predictions of clinical progression

main.delta = .1; % max death rate of tumor cells
main.gammaT = 1/6; % ec50 for phiD (on [0 .5]) inhibiting apoptosis
                      
main.daughters_reset = false; % if true, then daughter cells will have all monomers of FGFR3; otherwise will inherit concentrations from parent

main.extra = 15; % no longer need this so big because I will assume a neumann boundary condition; but bigger means more times using the same pde matrices for solving before the tumor grows beyond these bounds
              
%% FGFR3 signaling
main.V = main.cell_width^3; % m^3

main.RT = ((1e4/6.02214e23)*1e9/main.V)*1e-3; % in nM (per cell);
main.kf = 460.8; % in per nM per day
main.kr = 400; % dissociation rate in days^-1
main.kp = 112.32; % recycling rate of dimers in days^-1

%% anti-FGFR3 related
main.k_on_R = 1e2; % inhibitor binding to monomer (in per nM per day)
main.k_off_R = 7.8e2; % unbinding of inhibitor and monomer (in per day)
main.k_on_D = 1e2; % inhibitor binding to active dimer (in per nM per day)
main.k_off_D = 7.8e2; % unbinding of inhibitor and active dimer (in per day)

main.aFGFR3_influx = 10; % the rate that concentration is accumulated in periphery (in per day)
main.aFGFR3_eflux = 0; % rate that concentration leaves the periphery circulation (in per day)

main.aFGFR3_degradation = 144; % combo rate of anti-FGFR3 diffusing out of tumor + rate of degradation within tumor (per day)

main.inhibitor_entry = 'everywhere'; % where the inhibitor enters the TME

main.aFGFR3_circ0 = 1000; % initial concentration of circulating inhibitor in nM

main.aFGFR3_sysdecay = log(2)/.125; % decay rate of drug from system

%% pde related
main.aFGFR3_diffusion = 1e-5; % diffusion coefficient for fgfr3 inhibitor in TME (in m^2/d)
main.max_pde_dt = (.5/60) / 24;   % biggest allowable dt for updating the PDE
main.n_between = 0; % number of grid points between possible cell locations
main.deltaX = 1/(main.n_between+1); % number of grid points between possible cell locations in all directions; needs to be 1/n, for n\in\mathbb{Z}

%% vasculature (only for influx method 'multiple lines'
main.vessel_spacing = 200e-6*2; % all cells are within 200 micrometers of a blood vessel, so space them out 400 micrometers in x and y directions
%% construct struct
p = struct('plot_properties',plot_properties,'main',main);