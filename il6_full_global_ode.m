function dx = il6_full_global_ode(t,x,solver)

cell_bound_concentrations = reshape(x(1:solver.num_cell_bound),length(solver.bound_inds),solver.n_regions,solver.num_types);
free_concentrations = reshape(x(solver.num_cell_bound+1:end),[],solver.n_regions);
pk_concentrations = solver.pk_concentrations; % [circulation;periphery] concentrations for each pk substrate

dx_free = zeros(size(free_concentrations));
for si = 1:2
    conc_diffs = free_concentrations(si,:)'-free_concentrations(si,:);
    diffusion = 6*solver.substrate_pars(si).diffusion*sum(conc_diffs.*solver.M);
    dx_free(si,:) = dx_free(si,:) + diffusion - solver.substrate_pars(si).degradation*free_concentrations(si,:);
    if solver.is_pk(si)
        dx_free(si,:) = dx_free(si,:) + solver.regional_BV_prop .* solver.substrate_pars(si).fluid_exchange_rate .* (solver.substrate_pars(si).circulation_update_dual(t)*[0;pk_concentrations(:,si)]-free_concentrations(si,:));
    end
end

dx_bound = zeros(5,solver.n_regions,solver.num_types);

for ri = 1:solver.n_regions
    for ti = 1:solver.num_types
        ic = [cell_bound_concentrations(:,ri,ti);free_concentrations(:,ri)];
        dx_bound(:,ri,ti) = il6_ode([],ic(solver.perm_order),solver.agent_ode_pars);
    end
end

dx_free(1,:) = dx_free(1,:) + sum(dx_bound(2,:,:).*solver.type_proportions,3);
dx_free(2,:) = dx_free(2,:) + sum(dx_bound(4,:,:).*solver.type_proportions,3);

for si = 1:2
    if solver.is_dirichlet_at_blood_vessels(si)
        dx_free(si,solver.regional_BV_prop>0) = 0;
    end
end

dx = [reshape(dx_bound([1,3,5],:,:),[],1);dx_free(:)];