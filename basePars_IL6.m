function out = basePars_IL6()

flags.plotFigs = false;  % whether or not to plot any figures             
flags.makeMovie = false;   % whether or not to make a movie (need plotFigs=true)
flags.plot_every = 360; % how many minutes to wait before plotting again
flags.plot_offset = 0; % if plot_every>1, this determines which of the plot_every to plot
flags.plot_inhibitor = false; % if running the local method, will plot inhibitor after PDE step
flags.track_all_phiD = false; % whether or not to track all phiD values
flags.qs = linspace(0,1,11); % quantiles to compute for phiD
flags.go_without_tumor = false; % whether to continue even if the tumor is gone
flags.growing_me_size = false; % whether or not to allow the microenvironment grid to grow
flags.make3dplot = false; % whether or not to make a 3d plot during simulation

initialization_pars.initialize_cell_receptors_by_type = true; % whether the type of cell determines initial receptor concentrations
initialization_pars.start_with_drug = [false,false]; % whether or not to start with the drug in circulation
initialization_pars.circulation_concentration = [0,0];  % how much drug to start in circulation
initialization_pars.periphery_concentration = [0,0]; % how much drug to start in periphery

method = "local";

plot_properties.plotLocations = false; % whether or not to plot a 3D rendering of tumor (need plotFigs=true)
plot_properties.nrows = 4;
plot_properties.ncols = 4;
plot_properties.type_names = ["stem","progenitor","TD"];
inds.location_inds = 1:3; % locations of grid points
inds.subs_inds = 4:6; % subscripts of where tumor cells are in L array
inds.ind_ind = 7; % column containing the linear index of where cells are on the IL6 grid
inds.type_ind = 8; % column containing the type of tumor cell (1=stem cell, 2=progenitor, 3=terminally differentiated)
inds.event_ind = 9; % column containing the  event the cell is currently undergoing
inds.proliferation_timer_ind = 10; % column containing the minimum time to wait until next proliferation
inds.num_prolifs = 11; % tracking number of prolifs remaining for progenitor cells
inds.tumor_receptors_inds = 12:16; % columns containing the receptor information (unbound IL6R, free IL6 nearby, IL6R bound to IL6, free aIL6R nearby, IL6R bound to aIL6R)

simpars.censor_date = 60*24; % number of minutes to simulate

codensity_pars.track_intensities = false; % whether or not to save tracked.codensity_intensities (it can be a large variable)
codensity_pars.cod_edges = [1:1:22,Inf]; % quantiles to compute for phiD
codensity_pars.codensity_count = 21; % the nth nearest cell for codensity computations

pars.num_types = 3; % number of types of cells
initialization_pars.N0 = [10,5,0]; % initial tumor popoulations by type

pars.TME_size = [10,10,10]; % number of lattice points in grid

pars.cell_width = 20; % in microns; bladder cancer cell is about 20micrometers in diameter

pars.desired_dt = 30; % desired number of minutes per step

%% neighbor parameters
pars.occmax = 20; % below this threshold, a tumor cell can divide; at and above, too many neighbors and so doesn't proliferate
pars.neighbors = allCombos(-1:1,-1:1,-1:1,'matlab'); % local coordinates of neighbors in Moore neighborhood
pars.neighbors(all(pars.neighbors==0,2),:) = []; % don't count self as neighbor

%% tumor parameters
prolif_rate_S = 0.6/(24*60); % proliferation rate of cancer stem cells (per minute)
prolif_rate_E = (log(2)/1.04)/(24*60); % proliferation rate of progenitor cells (per minute)
prolif_rate_D = 0/(24*60); % proliferation rate of terminally differentiated cells (per minute)

event_pars.prolif_rate = [prolif_rate_S,prolif_rate_E,prolif_rate_D];

pars.min_prolif_wait = 9*60; % number of minutes all cells must wait at minimum between proliferations

update_pars.P_Smin_star = 0.014; % minimum probability of symmetric division for cancer stem cells (approached as number of cancer stem cells --> infty)
update_pars.P_Smax = 0.9; % max probability of symmetric division for cancer stem cells (approached as number of cancer stem cells --> 0)
update_pars.P_ec50 = 100; % number of stem cells at which the proliferation rate of stem cells is halfway between P_Smin_star and P_Smax
update_pars.hill_coeff = 2.6; % hill coefficient for suppression of stem cell self-renewal
update_pars.mu_S = 0.04; % modulation parameter for effect of IL6 on the minimum self-renewal rate of stem cells

% P: probability of stem cell self-renewal
% P_Smin: minimum probability of stem cell self-renewal
% S: number of stem cells
% phi_S: fractional occupancy of stem cells
% P_Smin = @(phi_S) mu_S * (P_Smax - P_Smin_star)*phi_S + P_Smin_star;
% P = @(S,phi_S) = (P_Smax - P_Smin(phi_S)) * P_ec50^hillcoeff / (P_ec50^hill_coeff+S^hill_coeff) + P_Smin_star;

delta_S = 1.5*0.6*0.014 / (24*60); % max death rate of cancer stem cells (per minute)
delta_E = 0.0612 / (24*60); % max death rate of progenitor cells (per minute)
delta_D = 0.0612 / (24*60); % max death rate of terminally differentiated cells (per minute)

event_pars.delta = [delta_S,delta_E,delta_D];

% ec50_S = 1/2.38; % ec50 for phiD inhibiting apoptosis for cancer stem cells
% ec50_E = 1/2.38; % ec50 for phiD inhibiting apoptosis for progenitor cells
% ec50_D = 1/2.38; % ec50 for phiD inhibiting apoptosis for terminally differentiated cells
% 
% event_pars.ec50 = [ec50_S,ec50_E,ec50_D];
event_pars.ec50 = 1/2.38;

update_pars.w = 2; % number of proliferations as progenitor cell before becoming terminally differentiated and no longer proliferating
event_pars.move_rate_in_microns = 2; % assume tumor cells move 2um/minute

update_pars.chemotaxis_substrate_ind = 1;
update_pars.chemotaxis_enabled = false(1,3);
update_pars.biased_move_probability = 0.5;

event_pars.event_prob_fn = @il6EventProbabilities;
%% IL6-IL6R signaling

RT_S = 1.66e-6; % femtomoles of IL6R on cancer stem cells
RT_E = 0.125*(1.66e-6); % femtomoles of IL6R on progenitor cells
RT_D = 0.125*(1.66e-6); % femtomoles of IL6R on terminally differentiated cells

pars.RT_nanomoles = [RT_S,RT_E,RT_D] * 1e-6; % number of nanomoles of IL6R on a tumor cell for all types
                           

pars.phiD_ind = 3; % index of receptors to compute phiD

substrate_pars.name = "IL6";
substrate_pars.diffusion = 1.08e2... % diffusion rate of IL12 (assume similar to IL6 even though it seems to be ~3x the molecular weight of IL6) in cm^2/day
                     *1e6... % convert cm^2 to um^2
                     /(24*60); % convert per day to per minute;
substrate_pars.degradation = 0.4152... % spontaneous degradation of IL6 in microenvironment (in per day)
              /(24*60); % in per minute
solver.agent_ode_pars.secretion_nanomoles = 7e-7... % rate of production of IL6 by each tumor cell (in femtomoles per day)
           *1e-6... % convert femtomoles to nanomoles
           /(24*60); % convert per day to per minute

substrate_pars.receptor_ind = 2; % IL6 is the second of the receptors on agents

solver.agent_ode_pars.kf = 2.72e4... % association rate of IL6 and IL6R (in per M per second)
          *1e-9*60; % in per nM per minute 
solver.agent_ode_pars.kr = 1.04e-2... % dissociation rate of IL6 and IL6R (in per second)
          *60; % in per minute
solver.agent_ode_pars.kp = 24.95... % internalization/recycling rate of IL6-IL6R complexes back to IL6R (in per day)
          /(24*60); % in per minute

substrate_pars(1).total_dual = [0;0;1;0;0]; % dual vector to take in receptor concentrations on agent and output total IL6 on agent
solver.is_dirichlet_at_blood_vessels(1) = true;
solver.dirichlet_condition(1) = 1; % assumed constant concentration of IL6 near blood vessels due to secretion by endothelial cells

%% anti-IL6R related

solver.agent_ode_pars.Kd_I = 2.54; % dissociation constant of aIL6R and IL6R (in nM)
solver.agent_ode_pars.kr_I = 21.6... % unbinding of inhibitor and IL6R (in per day)
            /(24*60); % in per minute

substrate_pars(2).name = "aIL6R";
substrate_pars(2).fluid_exchange_rate = 10; % the rate that the concentration differential between the blood and perivascular region decays (in per minute)
substrate_pars(2).k12 = 14.30... % rate that amount of aIL6R moves from circulation to periphery (in per day)
                /(24*60); % in per minute
substrate_pars(2).k21 = 5.55... % rate that amount of aIL6R moves from periphery to circulation (in per day)
                /(24*60);
substrate_pars(2).sysdecay = 0.004... % decay rate of aIL6R from circulation (in per day)
                     /(24*60); % in per minute
substrate_pars(2).receptor_ind = 4; % aIL6R is the fourth of the receptors on agents


pars.aIL6R_circ0 = 22; % initial concentration of circulating inhibitor in ug/mL


substrate_pars(2).degradation = 0; % rate of spontaneous degradation of aIL6R in microenvironment (in per minute)
substrate_pars(2).diffusion = 10.0... % diffusion coefficient for aIL6R in microenvironment (in um^2/second)
                      *60; % in um^2/minute

substrate_pars(2).total_dual = [0;0;0;0;1]; % dual vector to take in receptor concentrations on agent and output total aIL6R on agent
solver.is_dirichlet_at_blood_vessels(2) = false;

%% vasculature (only for influx method 'multiple lines')
pars.vessel_spacing = 200*2; % all cells are within 200 micrometers of a blood vessel, so space them out 400 micrometers in x and y directions

%% construct solver struct
solver.max_pde_dt = 0.5;   % biggest allowable dt for updating the PDE (in minutes)
solver.setup_done = false;
solver.num_substrates = 2;
solver.is_pk = [false,true]; % IL6 does not have PK dynamics, aIL6R does
solver.substrate_entry = "floor"; % where aIL6R enters the microenvironment
solver.plot_every_pde_step = false; % if running the local method, this can be used to plot every pde step to view microenvironment concentration
solver.is_present = [true,false]; % whether or not the substrates are present
solver.agent_ode = @il6_ode;
solver.global_method_full_ode = @il6_full_global_ode;
solver.substrate_pars = substrate_pars;
solver.separate_TOV_and_AMB = false; % whether to subdivide regions by being in TOV or AMB
solver.m = [1,1,1]; % size of regions to use in global method if substrate_entry is set to cell_source
%% construct struct
out = struct('method',method,...
    'plot_properties',plot_properties,...
    'pars',pars,...
    'flags',flags,...
    'simpars',simpars,...
    'initialization_pars',initialization_pars,...
    'solver',solver,...
    'inds',inds,...
    'event_pars',event_pars,...
    'update_pars',update_pars,...
    'codensity_pars',codensity_pars);

