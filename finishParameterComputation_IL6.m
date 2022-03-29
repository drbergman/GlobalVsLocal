function input = finishParameterComputation_IL6(input)

input.pars.cell_volume = input.pars.cell_width^3; % cell volume in micron^3 (assume cells are cubes)

input.event_pars.move_rate = input.event_pars.move_rate_in_microns/input.pars.cell_width; % convert to lattice edge lengths per minute

input.pars.RT = input.pars.RT_nanomoles / input.pars.cell_volume... % concentration of IL6R on a tumor cell for all types (nanomoles / um^3)
                           *1e15; % convert to nM (nanomoles/liter)

input.solver.agent_ode_pars.secretion = input.solver.agent_ode_pars.secretion_nanomoles/input.pars.cell_volume... % divide by volume of cell to get concentration in nanomoles/um^3/minute
           *1e15; % convert per um^3 to per L to get units of nM/minute

input.solver.agent_ode_pars.kf_I = input.solver.agent_ode_pars.kr_I/input.solver.agent_ode_pars.Kd_I; % inhibitor binding to IL6R (in per nM per minute)

input.solver.substrate_pars(2).circ_to_periph_volume_ratio = input.solver.substrate_pars(2).k21/input.solver.substrate_pars(2).k12; % assume that circulation and periphery concentrations are at quasi-equilibrium when their concentrations are equal


%% parameter computations
input.pars.Nsteps = ceil(input.simpars.censor_date/input.pars.desired_dt);
input.pars.dt = input.simpars.censor_date/input.pars.Nsteps; % actual time step

%% finish plot properties
input.plot_properties.num_types = input.pars.num_types;
input.plot_properties.censor_date = input.simpars.censor_date;

%% finish initialization pars
input.initialization_pars.growing_me_size = input.flags.growing_me_size;
input.initialization_pars.num_types = input.pars.num_types;
input.initialization_pars.TME_size = input.pars.TME_size;
input.initialization_pars.neighbors = input.pars.neighbors;
input.initialization_pars.prolif_rate = input.event_pars.prolif_rate;
input.initialization_pars.min_prolif_wait = input.pars.min_prolif_wait;
input.initialization_pars.w = input.update_pars.w; % number of proliferations as progenitor cell before becoming terminally differentiated and no longer proliferating
input.initialization_pars.initial_receptor_concentrations = [input.pars.RT',zeros(3,4)]; % columns: free IL6R, free IL6, IL6-IL6R, free aIL6R, aIL6R-IL6R


%% finish solver
input.solver.Tf = input.pars.dt;
input.solver.cell_width = input.pars.cell_width;
input.solver.method = input.method;
input.solver.num_types = input.pars.num_types;
input.solver.substrate_inds = [input.solver.substrate_pars.receptor_ind];
input.solver.bound_inds = setdiff(1:length(input.inds.tumor_receptors_inds),input.solver.substrate_inds);
input.solver.growing_me_size = input.flags.growing_me_size;

%% finish codensity_pars
if input.codensity_pars.track_intensities
    ne = length(input.codensity_pars.cod_edges);
    M = sparse(ne,ne-1);
    M(2:ne+1:end) = 1;
    M(1:ne+1:end) = 1;
    M = M./sum(M,2);
    input.codensity_pars.cod_intensity_mat = M;
end

%% create update_pars
input.update_pars.dt = input.pars.dt;
input.update_pars.occmax = input.pars.occmax;
input.update_pars.neighbors = input.pars.neighbors;
input.update_pars.min_prolif_wait = input.pars.min_prolif_wait;
input.update_pars.growing_me_size = input.flags.growing_me_size;
input.update_pars.RT = input.pars.RT;