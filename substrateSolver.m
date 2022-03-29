function [tumors,substrate,timer] = substrateSolver(tumors,substrate,timer,solver,inds)

if solver.method=="local"
    IC = tumors(:,inds.tumor_receptors_inds);
    NT = size(tumors,1);

    if solver.plot_every_pde_step
        figure(12345);
        ax = gca;
        plot(ax,1:solver.nz,[squeeze(min(substrate(1).me_concentration,[],1:2)),squeeze(max(substrate(1).me_concentration,[],1:2))])
        drawnow
    end

    for ti = 2:solver.nt
        for si = 1:solver.num_substrates
            if ~solver.is_present(si)
                continue;
            end
            if solver.is_pk(si)
                %% pk update

                substrate(si).me_concentration(solver.IE) = solver.substrate_pars(si).perivascular_update_dual*[reshape(substrate(si).me_concentration(solver.IE),1,[]);...
                    [substrate(si).circulation;substrate(si).periphery].*ones(1,sum(solver.IE,"all"))];
                pk_ic = [substrate(si).circulation;substrate(si).periphery];
                substrate(si).circulation = solver.substrate_pars(si).circulation_update_dual*[0;pk_ic];
                substrate(si).periphery = solver.substrate_pars(si).periphery_update_dual*[0;pk_ic];
            end

            %% diffusion decay terms
            timer.pde = tic;
            % x-diffusion
            for j = 1:solver.ny
                for k = 1:solver.nz
                    substrate(si).me_concentration(:,j,k) = solver.substrate_pars(si).Mx\substrate(si).me_concentration(:,j,k);
                end
            end

            % matlab solves this system faster if the one dimension varies along rows
            substrate(si).me_concentration = permute(substrate(si).me_concentration,[2,3,1]);

            % y-diffusion
            for i = 1:solver.nx
                for k = 1:solver.nz
                    %                 substrate(si).me_concentration(i,:,k) = solver.substrate_pars(si).My\reshape(substrate(si).me_concentration(i,:,k),[],1);
                    substrate(si).me_concentration(:,k,i) = solver.substrate_pars(si).My\substrate(si).me_concentration(:,k,i);
                end
            end

            % matlab solves this system faster if the one dimension varies along rows
            substrate(si).me_concentration = permute(substrate(si).me_concentration,[2,3,1]);

            % z-diffusion
            for i = 1:solver.nx
                for j = 1:solver.ny
                    %                 substrate(si).me_concentration(i,j,:) = solver.substrate_pars(si).Mz\reshape(substrate(si).me_concentration(i,j,:),[],1);
                    substrate(si).me_concentration(:,i,j) = solver.substrate_pars(si).Mz\substrate(si).me_concentration(:,i,j);
                end
            end

            % matlab solves this system faster if the one dimension varies
            % along rows (this step completes the 3-cycle and puts the
            % dimensions back in their original order)
            substrate(si).me_concentration = permute(substrate(si).me_concentration,[2,3,1]);

            if solver.is_dirichlet_at_blood_vessels(si)
                substrate(si).me_concentration(solver.IE) = solver.dirichlet_condition(si);
            end

            timer.pde_duration = timer.pde_duration + toc(timer.pde);
        end



        %% Cell supply/uptake
        for si = 1:solver.num_substrates
            IC(:,solver.substrate_pars(si).receptor_ind) = substrate(si).me_concentration(tumors(:,inds.ind_ind));
        end
        timer.ode = tic;

        for ci = NT:-1:1
            IC(ci,:) = iterativeEuler(IC(ci,:)',[],solver.dt,solver.agent_ode,solver.agent_ode_pars);
        end

        assert(all(IC>=0,"all"))

        if solver.plot_every_pde_step
            plot(ax,1:solver.nz,[squeeze(min(substrate(1).me_concentration,[],1:2)),squeeze(max(substrate(1).me_concentration,[],1:2))])
            drawnow
        end

        for si = 1:solver.num_substrates
            substrate(si).me_concentration(tumors(:,inds.ind_ind)) = IC(:,solver.substrate_pars(si).receptor_ind); % update locations that have a tumor cell

            if solver.is_dirichlet_at_blood_vessels(si)
                substrate(si).me_concentration(solver.IE) = solver.dirichlet_condition(si);
                if ti==solver.nt
                    IC(:,solver.substrate_pars(si).receptor_ind) = substrate(si).me_concentration(tumors(:,inds.ind_ind)); % reset this here at the last time step so tumor cells carry the correct free substrate concentrations
                end
            end
        end
        timer.ode_duration = timer.ode_duration+toc(timer.ode);

    end
    tumors(:,inds.tumor_receptors_inds) = IC;
elseif solver.method=="global"

    timer.ode = tic;

    tumors_in_region_by_type = cell(solver.n_regions,solver.num_types);
    IC = zeros(length(solver.bound_inds),solver.n_regions,solver.num_types);
    for ri = 1:solver.n_regions
        tumors_in_region = find(solver.regions(tumors(:,inds.ind_ind))==ri);
        for ti = 1:solver.num_types
            tumors_in_region_by_type{ri,ti} = tumors_in_region(tumors(tumors_in_region,inds.type_ind)==ti);
            if ~isempty(tumors_in_region_by_type{ri,ti})
                IC(:,ri,ti) = mean(tumors(tumors_in_region_by_type{ri,ti},inds.tumor_receptors_inds(solver.bound_inds)),1);
            else
                IC(:,ri,ti) = 0;
            end
            solver.type_proportions(1,ri,ti) = length(tumors_in_region_by_type{ri,ti})./solver.region_volumes(ri);
        end
    end
    substrate_IC = reshape([substrate.me_concentration],solver.n_regions,[])'; % [substrate.me_concentration] varies by first by region then by substrate, no matter the orientation of substrate(si).me_concentration
    IC = [IC(:);substrate_IC(:)]; % [substrate.me_concentration] varies by first by region then by substrate, no matter the orientation of substrate(si).me_concentration
    for si = 1:solver.num_substrates
        if solver.is_pk(si)
            solver.pk_concentrations(:,si) = [substrate(si).circulation;substrate(si).periphery];
        else
            solver.pk_concentrations(:,si) = [0;0];
        end
    end
    sol = ode15s(@(t,x) solver.global_method_full_ode(t,x,solver),[0 solver.dt],IC(:));
    cell_bound_concentrations = reshape(sol.y(1:solver.num_cell_bound,end),length(solver.bound_inds),solver.n_regions,solver.num_types);
    free_concentrations = reshape(sol.y(solver.num_cell_bound+1:end,end),[],solver.n_regions);
    for ri = 1:solver.n_regions
        for ti = 1:solver.num_types
            if ~isempty(tumors_in_region_by_type{ri,ti})
                tumors(tumors_in_region_by_type{ri,ti},inds.tumor_receptors_inds(solver.bound_inds)) = repmat(cell_bound_concentrations(:,ri,ti)',length(tumors_in_region_by_type{ri,ti}),1);
                tumors(tumors_in_region_by_type{ri,ti},inds.tumor_receptors_inds(solver.substrate_inds)) = repmat(free_concentrations(:,ri)',length(tumors_in_region_by_type{ri,ti}),1);
            end
        end
    end
    for si = 1:solver.num_substrates
        substrate(si).me_concentration = free_concentrations(si,:);
        if solver.is_pk(si)
            pk_ic = [substrate(si).circulation;substrate(si).periphery];
            substrate(si).circulation = solver.substrate_pars(si).circulation_update_dual(solver.dt)*[0;pk_ic];
            substrate(si).periphery = solver.substrate_pars(si).periphery_update_dual*[0;pk_ic];
        end
    end

    timer.ode_duration = timer.ode_duration+toc(timer.ode);

else
    error("unknown method")
end