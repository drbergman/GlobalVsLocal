function input = simTumor_IL6(input)

% TOV = tumor-occupied volume

Nsteps = input.pars.Nsteps;
go_without_tumor = input.flags.go_without_tumor;
if Nsteps == 0 || (isfield(input,"tumors") && isempty(input.tumors) && ~go_without_tumor) % nothing to simulate here
    return;
end

%% write down parameters (use structvars([structname]) to reproduce all this text
pars = input.pars;
flags = input.flags;
method = input.method;
solver = input.solver;
inds = input.inds;
event_pars = input.event_pars;
update_pars = input.update_pars;
codensity_pars = input.codensity_pars;
                    
track_all_phiD = flags.track_all_phiD;               

codensity_count = codensity_pars.codensity_count;

TME_size = pars.TME_size;
num_types = pars.num_types;
dt = pars.dt;
RT = pars.RT;

num_substrates = solver.num_substrates;

%% setup tumor and microenvironment
if ~isfield(input,"tumors")
    [tumors,grid,L] = initializeTumor_IL6(input.initialization_pars,...
        inds);
    solver.gridsize = grid.size;
    solver = setupSolver_IL6(solver);
    for si = num_substrates:-1:1
        if method=="local"
            substrate(si).me_concentration = zeros(grid.size);
        elseif method=="global"
            substrate(si).me_concentration = zeros(solver.n_regions,1);
        else
            error("unknown method")
        end
        if solver.is_pk(si)
            if input.initialization_pars.start_with_drug(si)
                substrate(si).circulation = input.initialization_pars.circulation_concentration(si);
                substrate(si).periphery = input.initialization_pars.periphery_concentration(si);
                solver.is_present(si) = true;
            else
                substrate(si).circulation = 0;
                substrate(si).periphery = 0;
            end
        end
        if solver.is_dirichlet_at_blood_vessels(si)
            if method=="local"
                substrate(si).me_concentration(solver.IE) = solver.dirichlet_condition(si);
                tumors(:,inds.tumor_receptors_inds(solver.substrate_inds(si))) = substrate(si).me_concentration(tumors(:,inds.ind_ind));
            elseif method=="global"
                substrate(si).me_concentration(solver.regional_BV_prop>0) = solver.dirichlet_condition(si);
                tumors(:,inds.tumor_receptors_inds(solver.substrate_inds(si))) = substrate(si).me_concentration(solver.regions(tumors(:,inds.ind_ind)));
            else
                error('unknown method')
            end
        end
    end
else
    tumors = input.tumors;
    grid = input.grid;
    L = input.L;
    substrate = input.substrate;
end

NT = size(tumors,1);


%% Tracking values
if isfield(input,'tracked')
    tracked = input.tracked;
    T = tracked.T(end);

    tum_apop = tracked.tum_apop(end,:);       % tracking total number of apoptosis by tumor cells
    
    names = fieldnames(tracked);
    counter = size(tracked.(names{1}),1)+1;
    for i = 1:length(names)
        if strcmp(names{i},'all_phiD')
            tracked.all_phiD = [tracked.all_phiD;cell(Nsteps,1)];
        else
            tracked.(names{i}) = cat(1,tracked.(names{i}),zeros([Nsteps,size(tracked.(names{i}),2:ndims(tracked.(names{i})))]));
        end
    end
else
    counter = 1;
    tracked.T    = zeros(Nsteps+1,1); % tracking the time
    tracked.NT   = zeros(Nsteps+1,num_types); % tracking tumor size by type
    tracked.phiD_mean = zeros(Nsteps+1,num_types); % mean phiD value for each type
    tracked.phiD_std = zeros(Nsteps+1,num_types); % std of phiD value for each type
    tracked.substrate_ambient = zeros(Nsteps+1,num_substrates); % tracking average ambient substrates in TME
    tracked.substrate_TOV = zeros(Nsteps+1,num_substrates); % tracking average TOV substrates in TME
    tracked.substrate_circ = NaN(Nsteps+1,num_substrates); % tracking concentration of substrates in circulation
    tracked.substrate_periph = NaN(Nsteps+1,num_substrates); % tracking concentration of substrates in periphery
    tracked.tum_apop = zeros(Nsteps+1,num_types);
    tracked.substrate_tum = zeros(Nsteps+1,num_types,num_substrates); % average (total) concentration of substrates bound to tumor cells
    tracked.mean_z = zeros(Nsteps+1,num_types); % average z coordinate of each type (useful when blood vessels on floor
    tracked.std_z = zeros(Nsteps+1,num_types); % average z coordinate of each type (useful when blood vessels on floor

    tracked.codensity = zeros(Nsteps+1,num_types); % average codensity by cell type
    tracked.codensity_std = zeros(Nsteps+1,num_types); % std of average codensity by cell type
    
    tracked.codensity_by_type = zeros(Nsteps+1,num_types,num_types); % average codensity with respect to each type by cell type

    tracked.substrate_by_z = zeros(Nsteps+1,TME_size(3),num_substrates);
    tracked.substrate_by_z_std = zeros(Nsteps+1,TME_size(3),num_substrates);

    tracked.wall_time = zeros(Nsteps,1);
    tracked.wall_time_ode = zeros(Nsteps,1);
    tracked.wall_time_pde = zeros(Nsteps,1);
    
    if track_all_phiD
        tracked.all_phiD = cell(Nsteps,num_types); % tracking all phi_d values for each cell
    end

    T = 0; % tracking time
    tum_apop = zeros(1,num_types);    % tracking total number of apoptosis of tumor cells by type
    
    %% set values at t=0
    tracked.NT(1,:) = input.initialization_pars.N0;
    for ti = 1:num_types % inds.type_ind
        phid_temp = tumors(tumors(:,inds.type_ind)==ti,inds.tumor_receptors_inds(3))/RT(ti);
        tracked.phiD_mean(1,ti) = mean(phid_temp);
        tracked.phiD_std(1,ti) = std(phid_temp);
        z_vals = (tumors(tumors(:,inds.type_ind)==ti,inds.subs_inds(3))-1)/grid.size(3);
        tracked.mean_z(1,ti) = mean(z_vals);
        tracked.std_z(1,ti) = std(z_vals);
    end

    [tracked.codensity(1,:),tracked.codensity_std(1,:),...
        tracked.codensity_by_type(1,:,:)] = ...
        codensityFunction(tumors(:,inds.location_inds),codensity_count,tumors(:,inds.type_ind),num_types);

    tracked = saveSubstrateInfo(tracked,solver,method,substrate,tumors(:,inds.ind_ind),counter);
    
    
    tracked.tum_apop(1,:) = tum_apop;
    counter = 2;
end

%% iterations
for i = 1:Nsteps
    timer.update = tic;
    
    timer.ode_duration = 0;
    timer.pde_duration = 0;
    
    if method=="local"

        %% solve PDE

        [tumors,substrate,timer] = substrateSolver_IL6(tumors,substrate,timer,solver,inds);

    elseif method=="global"
        [tumors,substrate,timer] = substrateSolver_IL6(tumors,substrate,timer,solver,inds);
    else
        error('unknown method')
    end % end of how to update tumors before cell fate decisions
    
    phiD = tumors(:,inds.tumor_receptors_inds(pars.phiD_ind))./RT(tumors(:,inds.type_ind))';
    for ti = num_types:-1:1 % inds.type_ind
        phiD_temp{ti} = phiD(tumors(:,inds.type_ind)==ti);
        tracked.phiD_mean(counter,ti) = mean(phiD_temp{ti});
        tracked.phiD_std(counter,ti) = std(phiD_temp{ti});
    end
    
    if track_all_phiD 
        tracked.all_phiD{counter-1,:} = phiD_temp;
    end
    
    p_matrix_tum = event_pars.event_prob_fn(tumors,event_pars,dt,inds,phiD);

    tumors(:,inds.event_ind) = arrayfun(@(ti) find(rand()<cumsum([p_matrix_tum(ti,:),1]),1),1:NT);
    
    for ti = 1:num_types % type index
        tum_apop(ti) = tum_apop(ti) + sum(tumors(tumors(:,inds.type_ind)==ti,inds.event_ind)==2);
    end
    
    active_ind = find(tumors(:,inds.event_ind)<=size(p_matrix_tum,2));
    order = active_ind(randperm(length(active_ind)));
    
    if method=="local"
        [tumors,L] = updateTumor_IL6(tumors,L,substrate,order,inds,grid,update_pars,method,phiD,solver);
            
    elseif method=="global"
        [tumors,L] = updateTumor_IL6(tumors,L,substrate,order,inds,grid,update_pars,method,phiD,solver);
    else
        error('unknown method')
    end
               
    %% clean up tumor stuff
    apoptosis_log = tumors(:,inds.event_ind)==2;
    L(tumors(apoptosis_log,inds.ind_ind)) = false;
    
    tumors(apoptosis_log,:)=[]; % get rid of dead cells
    tumors(tumors(:,inds.event_ind)>1,inds.proliferation_timer_ind) = max(0,tumors(tumors(:,inds.event_ind)>1,inds.proliferation_timer_ind)-dt); % update time until tumor cells can proliferate for all that did not proliferate (1) or were just created (0) because those timers were updated in updateTumor
    
    NT = size(tumors,1);
    
    T = T+dt;
    
    %% Update remaining tracked values
    tracked.T(counter,1)    = T;
    tracked.NT(counter,:)   = accumarray(tumors(:,inds.type_ind),1,[3,1]);
    tracked = saveSubstrateInfo(tracked,solver,method,substrate,tumors(:,inds.ind_ind),counter);
    tracked.tum_apop(counter,:) = tum_apop;
    for si = 1:num_substrates
        tracked.substrate_tum(counter,:,si) = accumarray(tumors(:,inds.type_ind),tumors(:,inds.tumor_receptors_inds)*solver.substrate_pars(si).total_dual,[3,1])'./tracked.NT(counter,:);
    end

    for ti = 1:num_types
        z_vals = (tumors(tumors(:,inds.type_ind)==ti,inds.subs_inds(3))-1)/grid.size(3);
        tracked.mean_z(counter,ti) = mean(z_vals);
        tracked.std_z(counter,ti) = std(z_vals);
    end

    [tracked.codensity(counter,:),tracked.codensity_std(counter,:),...
        tracked.codensity_by_type(counter,:,:)] = ...
        codensityFunction(tumors(:,inds.location_inds),codensity_count,tumors(:,inds.type_ind),num_types);

    tracked.wall_time_ode(counter-1,1) = timer.ode_duration;
    tracked.wall_time_pde(counter-1,1) = timer.pde_duration;
    tracked.wall_time(counter-1,1) = toc(timer.update);
    
    %% increase counter
    counter = counter+1;
    
    %% tumor gone?
    if NT == 0 && ~go_without_tumor % if all the tumor cells die, break the loop
        break
    end
    
end %%end of for

input.pars = pars;
input.flags = flags;
input.method = method;
input.solver = solver;
input.inds = inds;
input.event_pars = event_pars;
input.update_pars = update_pars;
input.grid = grid;
input.substrate = substrate;
input.tracked = tracked;
input.tumors = tumors;
input.L = L;
