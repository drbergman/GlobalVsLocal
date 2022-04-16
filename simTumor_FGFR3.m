function [tracked,tumors,aFGFR3_tracker,grids] = ...
    simTumor_FGFR3(tumors,pars,Nsteps,aFGFR3_tracker,grids,tracked)

% TOV = tumor-occupied volume

if Nsteps == 0 || (isempty(tumors) && ~pars.main.waitForImmune) % nothing to simulate here
    return;
end

location_inds = 1:3;
event_ind = 5;
proliferation_timer_ind = 6;
subs_inds = 7:9;
ind_ind = 10;

tumor_clearance_inds = 11:12;
tumor_receptors_inds = 14:19; % column indices in which the receptor information is stored

%% write down parameters
main = pars.main;

track_all_phiD = main.track_all_phiD;
dt = main.dt;
cell_width = main.cell_width;
alpha1 = main.alpha1;
alpha2 = main.alpha2;
delta = main.delta;
gammaT = main.gammaT;
extra = main.extra;
RT = main.RT;
aFGFR3_influx = main.aFGFR3_influx;
aFGFR3_eflux = main.aFGFR3_eflux;
aFGFR3_degradation = main.aFGFR3_degradation;
aFGFR3_sysdecay = main.aFGFR3_sysdecay;
inhibitor_entry = main.inhibitor_entry;
method = main.method;
aFGFR3_diffusion = main.aFGFR3_diffusion;
max_pde_dt = main.max_pde_dt;
deltaX = main.deltaX;
all_receptors_without_inhibitor = main.all_receptors_without_inhibitor;
kf = main.kf;                                                               
kr = main.kr;                                                               
kp = main.kp;                                                               
k_on_R = main.k_on_R;                                                       
k_off_R = main.k_off_R;                                                     
k_on_D = main.k_on_D;                                                       
k_off_D = main.k_off_D; 
vessel_spacing = main.vessel_spacing;

%% parameter-dependent setup begins here
if numel(tumors)==1
    %% preparing the tumor
    tumors = initializeTumor_FGFR3(tumors,RT,...
        location_inds,event_ind,tumor_clearance_inds,...
        tumor_receptors_inds,proliferation_timer_ind);
end

NT = size(tumors,1);

%% Set up grid
mins = min(tumors(:,location_inds),[],1);
maxs = max(tumors(:,location_inds),[],1);

if isempty(grids)
    xx = (mins(1)-extra):deltaX:(maxs(1)+extra); % I can make extra depend on the size of the tumor and that dimension
    yy = (mins(2)-extra):deltaX:(maxs(2)+extra);
    zz = (mins(3)-extra):deltaX:(maxs(3)+extra);    
else
    xx = grids{1};
    yy = grids{2};
    zz = grids{3};
end
grid_size = [length(xx),length(yy),length(zz)];
V_tot = prod(grid_size);

% (tumor location relative to start of lattice)/(length of each step) =
% num of lattice "right" of start; add 1 because start has index 1
tumors(:,subs_inds) = (tumors(:,location_inds)-repmat([xx(1),yy(1),zz(1)],[NT,1]))/deltaX + 1;
tumors(:,ind_ind) = sub2ind(grid_size,tumors(:,subs_inds(1)),tumors(:,subs_inds(2)),tumors(:,subs_inds(3))); % linear indices

aFGFR3_circ = aFGFR3_tracker{1};
aFGFR3 = aFGFR3_tracker{2};

if aFGFR3_circ>0 || any(aFGFR3>0)% set up fgfr3 inhibitor stuff
    if strcmp(method,'Local')
        aFGFR3_pde_constructs = struct('aFGFR3_degradation',aFGFR3_degradation,...
            'aFGFR3_sysdecay',aFGFR3_sysdecay,...
            'aFGFR3_diffusion',aFGFR3_diffusion,'t_pde',dt,...
            'nt',max(2,ceil(dt/max_pde_dt)),'aFGFR3_influx',aFGFR3_influx,...
            'aFGFR3_eflux',aFGFR3_eflux,...
            'inhibitor_entry',inhibitor_entry,'setup_done',false);
                        
        if numel(aFGFR3)<prod(grid_size)
            aFGFR3=aFGFR3+zeros(length(xx),length(yy),length(zz));
        end
    end
end

%% Tracking values
if ~isempty(tracked)
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
    counter = 2;
    tracked.T    = zeros(Nsteps+1,1); % tracking the time
    tracked.NT   = zeros(Nsteps+1,1); % tracking tumor size
    tracked.phiD_mean = zeros(Nsteps+1,1); % mean phiD value
    tracked.phiD_std = zeros(Nsteps+1,1); % std of phiD value
    tracked.aFGFR3_ambient = zeros(Nsteps+1,1); % tracking ambient aFGFR3 in TME
    tracked.aFGFR3_TOV = zeros(Nsteps+1,1);
    tracked.aFGFR3_circ = zeros(Nsteps+1,1);
    tracked.tum_apop = zeros(Nsteps+1,1);
    tracked.aFGFR3_tum = zeros(Nsteps+1,1); % average concentration of inhibitor on tumor cells
    
    tracked.expectedGrowth = zeros(Nsteps,1);
    tracked.simple_expectedGrowth = zeros(Nsteps,1);
    
    tracked.wall_time = zeros(Nsteps,1);
    tracked.wall_time_ode = zeros(Nsteps,1);
    tracked.wall_time_pde = zeros(Nsteps,1);
    
    if track_all_phiD
        tracked.all_phiD = cell(Nsteps,1);
    end
    T = 0; % tracking time
    
    tum_apop = 0;    % tracking total number of apoptosis by tumor cells
    
    %% set values at t=0
    tracked.NT(1) = NT;
    
    tracked.phiD_mean(1,:) = mean(tumors(:,tumor_receptors_inds(2))/RT);
    tracked.phiD_std(1,:) = std(tumors(:,tumor_receptors_inds(2))/RT);
    
    if any(aFGFR3>0,'all') || aFGFR3_circ>0
        tracked.aFGFR3_circ(1,1) = aFGFR3_circ;
        if strcmp(method,'Local')
            aFGFR3_ambient_temp = aFGFR3;
            aFGFR3_ambient_temp(tumors(:,ind_ind)) = [];
            
            tracked.aFGFR3_ambient(1,1) = mean(aFGFR3_ambient_temp);
            tracked.aFGFR3_TOV(1,1) = mean(aFGFR3(tumors(:,ind_ind)));
        else
            tracked.aFGFR3_ambient(1,1) = aFGFR3;
            tracked.aFGFR3_TOV(1,1) = mean(tumors(tumors(:,event_ind)>1,tumor_receptors_inds(3)));
        end
    end
    
    tracked.tum_apop(1,:) = tum_apop;
        
end

%% iterations
for i = 1:Nsteps
    timer.update = tic;
    
    timer.ode_duration = 0;
    timer.pde_duration = 0;
    
    if aFGFR3_circ==0 && all(aFGFR3==0,'all') % then new cells can be updated predictably and there's no PDE to solve
        tumors(:,tumor_receptors_inds) = repmat(all_receptors_without_inhibitor,NT,1);
    elseif strcmp(method,'Local') % each cell will be different and so will solve the ODE for 1/2 the interval, update the pde, then solve ODE for second half
        %% solve ODE once
        tumors(:,tumor_receptors_inds(3)) = aFGFR3(tumors(:,ind_ind));
        
        %% solve PDE
        [aFGFR3,aFGFR3_circ,aFGFR3_pde_constructs,tumors,timer] = ...
            substrateSolver_FGFR3(xx*cell_width,yy*cell_width,zz*cell_width,...
            aFGFR3,aFGFR3_circ,aFGFR3_pde_constructs,tumors,tumor_receptors_inds,ind_ind,timer,...
            kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D,inhibitor_entry,location_inds,cell_width,vessel_spacing);
        
        if pars.plot_properties.plot_inhibitor
            figure(54321);
            [sx,sy,sz,sc] = sphere_with_carveout(xx,yy,zz,aFGFR3);
            surf(sx,sy,sz,sc,'FaceColor','interp','EdgeColor','interp');
            title(sprintf('Time = %3.2f hours',T*24))
            view(144,40)
            drawnow
        end

    else 
        %% global method for update
        main.volume_prop = NT/V_tot;
        main.prop_tum_neighbors_that_are_tumors = main.volume_prop;
        main.prop_free_neighbors_that_are_free = 1 - main.volume_prop;

        switch inhibitor_entry
            case 'everywhere'
                main.prop_tov_receiving_inhibitor = 1;
                main.prop_amb_receiving_inhibitor = 1;
            case 'outside' % all locations outside tumor will receive inhibitor
                com = mean(tumors(:,location_inds),1); % tumor center of mass
                max_distance = max(sum((tumors(:,location_inds)-com).^2,2));
                main.prop_tov_receiving_inhibitor = 0;
                main.prop_amb_receiving_inhibitor = sum(sum((allCombos(xx,yy,zz)-com).^2,2)>max_distance,'all')/(V_tot-NT);
            case 'border'
                tumors_at_border = sum(any((tumors(:,location_inds)==(mins-extra)) | (tumors(:,location_inds)==(maxs+extra)),2));
                main.prop_tov_receiving_inhibitor = tumors_at_border/NT;
                main.prop_amb_receiving_inhibitor = (V_tot - prod(grid_size-2)-tumors_at_border)/(V_tot-NT);
            case 'line' % the line will be through x=y=0
                tumors_on_line = sum(all(tumors(:,location_inds(1:2))==0,2));
                main.prop_tov_receiving_inhibitor = tumors_on_line/NT;
                main.prop_amb_receiving_inhibitor = (length(zz)-tumors_on_line)/(V_tot-NT); 
            case 'plane' % the plane will be through x=0
                tumors_on_plane = sum(all(tumors(:,location_inds(1))==0,2));
                main.prop_tov_receiving_inhibitor = tumors_on_plane/NT;
                main.prop_amb_receiving_inhibitor = (length(yy)*length(zz)-tumors_on_plane)/(V_tot-NT);
            case 'multiple lines' % lines will be in z direction spaced out on a lattice
                xloc_bv = vessel_spacing*(ceil(xx(1)*cell_width/vessel_spacing):floor(xx(end)*cell_width/vessel_spacing));
                yloc_bv = vessel_spacing*(ceil(yy(1)*cell_width/vessel_spacing):floor(yy(end)*cell_width/vessel_spacing));
                [~,xind_bv] = min(abs(cell_width*xx-xloc_bv'),[],2);
                [~,yind_bv] = min(abs(cell_width*yy-yloc_bv'),[],2);
                tumors_on_a_line = sum(any(tumors(:,location_inds(1))==xx(xind_bv),2) & any(tumors(:,location_inds(2))==yy(yind_bv),2));
                main.prop_tov_receiving_inhibitor = tumors_on_a_line/NT;
                main.prop_amb_receiving_inhibitor = (length(zz)*length(xloc_bv)*length(yloc_bv)-tumors_on_a_line)/(V_tot-NT);
        end
        if NT > 0
            IC = [NT,0,mean(tumors(:,tumor_receptors_inds),1),aFGFR3]; % use average concentrations around each cell
            
            timer.ode = tic;
            [new_concentrations,aFGFR3] = fullGlobalODE_FGFR3(IC,main,dt,aFGFR3_circ);
            timer.ode_duration = toc(timer.ode);
            
            aFGFR3_circ = aFGFR3_circ * exp(-aFGFR3_sysdecay*dt);
            tumors(:,tumor_receptors_inds) = repmat(new_concentrations,[NT,1]);
        end
    end % end of how to update tumors before cell fate decisions
    
    phiD = tumors(:,tumor_receptors_inds(2))/RT;
    
    tracked.phiD_mean(counter,:) = mean(phiD);
    tracked.phiD_std(counter,:) = std(phiD);
    tracked.simple_expectedGrowth(counter-1) = mean(alpha1+alpha2*phiD - delta./(1+phiD/gammaT));
    
    if track_all_phiD 
        tracked.all_phiD{counter-1} = phiD; % do this here before it's changed later
    end
        
    p_matrix_tum = dt*[max(0,dt-tumors(:,proliferation_timer_ind))/dt.*(alpha1+alpha2*phiD),... % proliferation probabilities; want to make this density-dependent; any with timer>dt must wait, any with timer <=dt has their probability decreased depending on this timer
        delta./(1+phiD/gammaT)]; % spontaneous apoptosis probabilities
    
    tracked.expectedGrowth(counter-1) = mean(p_matrix_tum,1)*[1;-1];
    
    assert(all(sum(p_matrix_tum,2)<=1)) % make sure the probabilities do not add up to more than 1
    
    tumors(:,event_ind) = arrayfun(@(ti) find(rand()<cumsum([p_matrix_tum(ti,:),1]),1),1:NT);
    
    tum_apop = tum_apop + sum(tumors(:,event_ind)==2);
    
    active_ind = find(tumors(:,event_ind)<=size(p_matrix_tum,2));
    order = active_ind(randperm(length(active_ind)));
    
    switch method
        case 'Local' % the local method does not update aFGFR3 here and has no need for V_tot
            [tumors,~]      = updateTumor_FGFR3(main,order,tumors,...
                1    ,NT,0,location_inds,... % set aFGFR3 to 0 because will set new inhibitor concentration right before solving ODE first time
                event_ind,proliferation_timer_ind,tumor_receptors_inds);
            
        case 'Global'
            [tumors,aFGFR3] = updateTumor_FGFR3(main,order,tumors,...
                V_tot,NT,aFGFR3,location_inds,...
                event_ind,proliferation_timer_ind,tumor_receptors_inds);
    end
    
    %% update locations on grid
    tumors(NT+1:end,subs_inds) = (tumors(NT+1:end,location_inds)-repmat([xx(1),yy(1),zz(1)],[size(tumors,1)-NT,1]))/deltaX + 1;
    tumors(NT+1:end,ind_ind) = sub2ind(grid_size,tumors(NT+1:end,subs_inds(1)),tumors(NT+1:end,subs_inds(2)),tumors(NT+1:end,subs_inds(3))); % linear indices
    
    if strcmp(method,'Local') && ndims(aFGFR3)==3 % assign inhibitor concentration to new tumor cells
        tumors(NT+1:end,tumor_receptors_inds(3)) = aFGFR3(tumors(NT+1:end,ind_ind));
    end
        
    %% clean up tumor stuff
    apoptosis_log = tumors(:,event_ind)==2;
    
    tumors(apoptosis_log,:)=[]; % get rid of dead cells
    tumors(tumors(:,event_ind)>1,proliferation_timer_ind) = max(0,tumors(tumors(:,event_ind)>1,proliferation_timer_ind)-dt); % update time until tumor cells can proliferate for all that did not proliferate (1) or were just created (0) because those timers were updated in updateTumor
    
    NT = size(tumors,1);
    
    mins = min(tumors(:,location_inds),[],1);
    maxs = max(tumors(:,location_inds),[],1);
    
    T = T+dt;
    
    %% Update remaining tracked values
    tracked.T(counter,1)    = T;
    tracked.NT(counter,1)   = NT;
    if any(aFGFR3>0,'all') || aFGFR3_circ>0
        tracked.aFGFR3_circ(counter,1) = aFGFR3_circ;
        if strcmp(method,'Local')
            aFGFR3_ambient_temp = aFGFR3;
            aFGFR3_ambient_temp(tumors(:,ind_ind)) = [];
            
            tracked.aFGFR3_ambient(counter,1) = mean(aFGFR3_ambient_temp);
            tracked.aFGFR3_TOV(counter,1) = mean(aFGFR3(tumors(:,ind_ind)));
        else
            tracked.aFGFR3_ambient(counter,1) = aFGFR3;
            tracked.aFGFR3_TOV(counter,1) = mean(tumors(:,tumor_receptors_inds(3)));
        end
    end
    tracked.tum_apop(counter,:) = tum_apop;
    tracked.aFGFR3_tum(counter,:) = mean(tumors(:,tumor_receptors_inds)*[0;0;0;1;2;1]);
    
    tracked.wall_time_ode(counter-1,1) = timer.ode_duration;
    tracked.wall_time_pde(counter-1,1) = timer.pde_duration;
    tracked.wall_time(counter-1,1) = toc(timer.update);
    
    %% increase counter
    counter = counter+1;
    
    %% tumor gone?
    if NT == 0 % if all the tumor cells die, break the loop
        break
    end
    
end %%end of for

% save these for this patient to be used in subsequent calls to
% simTumor_FGFR3
aFGFR3_tracker = {aFGFR3_circ,aFGFR3};

grids = {xx,yy,zz};