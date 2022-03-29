function [tumors,grid,L] = initializeTumor_IL6(pars,inds)

ind_names = fieldnames(inds);
NT = sum(pars.N0);
NT_cum = [0,cumsum(pars.N0)];
max_ind = 1;
for fni = 1:length(ind_names)
    max_ind = max([max_ind,inds.(ind_names{fni})]);
end

tumors = zeros(NT,max_ind);


tum_inds = randsample(1:prod(pars.TME_size),NT_cum(end),false);
if isfield(inds,'ind_ind')
    tumors(:,inds.ind_ind) = tum_inds;
end
if isfield(inds,'subs_inds') || isfield(inds,'location_inds')
    [tum_subs(:,1),tum_subs(:,2),tum_subs(:,3)] = ind2sub(pars.TME_size,tum_inds);
    if isfield(inds,'subs_inds')
        tumors(:,inds.subs_inds) = tum_subs;
    end

    if isfield(inds,'location_inds')
        tumors(:,inds.location_inds) = tum_subs-floor(0.5*(pars.TME_size+1)); % move center of TME to (0,0,0) (or at least near that) reason: if blood vessels are created based on location, this will help ensure uniformity across tumor sizes and random starts to the tumor
    end
end

grid.size = pars.TME_size;
grid.xx = (1:grid.size(1))-floor(0.5*(pars.TME_size(1)+1));
grid.yy = (1:grid.size(2))-floor(0.5*(pars.TME_size(2)+1));
grid.zz = (1:grid.size(3))-floor(0.5*(pars.TME_size(3)+1));

grid.V_tot = prod(grid.size);
grid.rel_pos_ind = pars.neighbors*cumprod([1,grid.size(1:2)])'; % for Moore neighborhood


for ti = 1:pars.num_types
    tumors(NT_cum(ti)+1:NT_cum(ti+1),inds.type_ind) = ti;
end
tumors(NT_cum(2)+1:NT_cum(3),inds.num_prolifs) = pars.w;


tumors(:,inds.event_ind) = 0; % event index

if pars.initialize_cell_receptors_by_type
    tumors(:,inds.tumor_receptors_inds) = pars.initial_receptor_concentrations(tumors(:,inds.type_ind),:);
else
    tumors(:,inds.tumor_receptors_inds) = pars.initial_receptor_concentrations;
end

%% determining remaining wait times for proliferation on tumor cells (see my Onenote document in ABM>Rates--Exponential>Going backwards in time)
if length(pars.prolif_rate)>1
    prolif_rate = pars.prolif_rate(tumors(:,inds.type_ind));
    prolif_rate = reshape(prolif_rate,[],1); % ensure is column vector
else
    prolif_rate = pars.prolif_rate;
end
u = rand(NT,1);
i1 = u>=1./(1+pars.min_prolif_wait*prolif_rate);
i2 = ~i1;
x = zeros(NT,1);
x(i1) = u(i1).*(1+pars.min_prolif_wait*prolif_rate(i1))./prolif_rate(i1) - 1./prolif_rate(i1) - pars.min_prolif_wait;
x(i2) = log((1+pars.min_prolif_wait*prolif_rate(i2)).*u(i2))./prolif_rate(i2) - pars.min_prolif_wait;
tumors(:,inds.proliferation_timer_ind) = max(0,x+pars.min_prolif_wait); % time until next possible proliferation (days)

assert(all(tumors(:,inds.subs_inds(3))>=1))

%% initialize lattice
L = false(grid.size);
L(tumors(:,inds.ind_ind)) = true;