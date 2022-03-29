function solver = setupSolver(solver)

if solver.method=="global"
    switch solver.substrate_entry
        case "floor"
            % make regions based on height from floor
            solver.regions = reshape(repelem(1:solver.gridsize(3),prod(solver.gridsize(1:2))),solver.gridsize);
            solver.n_regions = solver.gridsize(3);
        case "cell_source"
            solver.regions = zeros(solver.gridsize);
            region_dims = solver.m;
            region_grid_size = solver.gridsize./region_dims;
            solver.regions = repelem((1:region_grid_size(1))',solver.m(1)) + ...
                repelem((0:region_grid_size(2)-1)*region_grid_size(1),solver.m(2)) + ...
                reshape(repelem(0:region_grid_size(3)-1,solver.m(3)),1,1,[])*prod(region_grid_size(1:2));
            solver.n_regions = prod(solver.gridsize./solver.m);
        otherwise
            warning('this inhibitor entry not yet completed, defaulting to one region')
            solver.regions = ones(solver.gridsize);
            solver.n_regions = 1;
    end
end

if any(solver.is_pk)
    switch solver.substrate_entry
        case "everywhere"
            solver.IE = true(solver.gridsize);
        case "border"
            solver.IE = true(solver.gridsize);
            solver.IE(2:end-1,2:end-1,2:end-1) = false;
        case "line"
            if solver.growing_me_size
                warning("If the grid is expanding, the location of the blood vessels may change.")
            end
            solver.IE = false(solver.gridsize);
            xind = round(0.5*solver.gridsize(1));
            yind = round(0.5*solver.gridsize(2));
            solver.IE(xind,yind,:) = true;
        case "plane"
            if solver.growing_me_size
                warning("If the grid is expanding, the location of the blood vessels may change.")
            end
            solver.IE = false(solver.gridsize);
            xind = round(0.5*solver.gridsize(1));
            solver.IE(xind,:,:) = true;
        case "floor"
            if solver.growing_me_size
                warning("If the grid is expanding, the location of the blood vessels may change.")
            end
            solver.IE = false(solver.gridsize);
            solver.IE(:,:,1) = true;
        case "outside"
            error("blood vessels outside the tumor not implemented")
        case "multiple lines"
            warning("If the grid is expanding, the location of the blood vessels may change.")
            solver.IE = false(solver.gridsize);
            x_center_ind = round(0.5*solver.gridsize(1));
            y_center_ind = round(0.5*solver.gridsize(2));
            vessel_ind_x = round(unique([x_center_ind:(-solver.vessel_spacing/solver.cell_width):0,x_center_ind:(solver.vessel_spacing/solver.cell_width):solver.gridsize(1)+1]));
            vessel_ind_y = round(unique([y_center_ind:(-solver.vessel_spacing/solver.cell_width):0,y_center_ind:(solver.vessel_spacing/solver.cell_width):solver.gridsize(2)+1]));
            vessel_ind_x(vessel_ind_x<1 | vessel_ind_x>solver.gridsize(1)) = [];
            vessel_ind_y(vessel_ind_y<1 | vessel_ind_y>solver.gridsize(2)) = [];
            solver.IE(vessel_ind_x,vessel_ind_y,:) = true;
        case "cell_source"
            solver.IE = false(solver.gridsize);
        otherwise
            error("Have not added ""%s"" to the local solver.method\n",solver.substrate_entry)
    end
end




if solver.method=="local"
    solver.nt = ceil(solver.Tf/solver.max_pde_dt)+1;
    solver.dt = solver.Tf/solver.nt;

    solver.nx = solver.gridsize(1);
    solver.ny = solver.gridsize(2);
    solver.nz = solver.gridsize(3);

    for si = 1:solver.num_substrates
        %% set up in x direction
        r_ind_ux = 1:solver.nx-1;
        c_ind_ux = 2:solver.nx;
        v_ux = -solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2 * ones(1,solver.nx-1);

        r_ind_dx = 1:solver.nx;
        c_ind_dx = 1:solver.nx;
        v_dx(1,2:solver.nx-1) = 1 + solver.dt*solver.substrate_pars(si).degradation/3 + 2*solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2;
        v_dx([1,solver.nx]) = 1 + solver.dt*solver.substrate_pars(si).degradation/3 + solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2;

        r_ind_lx = 2:solver.nx;
        c_ind_lx = 1:solver.nx-1;
        v_lx = -solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2 * ones(1,solver.nx-1);

        solver.substrate_pars(si).Mx = sparse([r_ind_ux,r_ind_dx,r_ind_lx],...
            [c_ind_ux,c_ind_dx,c_ind_lx],...
            [v_ux    ,v_dx   ,v_lx     ],solver.nx,solver.nx);

        %% set up in y direction
        r_ind_uy = 1:solver.ny-1;
        c_ind_uy = 2:solver.ny;
        v_uy = -solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2 * ones(1,solver.ny-1);

        r_ind_dy = 1:solver.ny;
        c_ind_dy = 1:solver.ny;
        v_dy(1,2:solver.ny-1) = 1 + solver.dt*solver.substrate_pars(si).degradation/3 + 2*solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2;
        v_dy([1,solver.ny]) = 1 + solver.dt*solver.substrate_pars(si).degradation/3 + solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2;

        r_ind_ly = 2:solver.ny;
        c_ind_ly = 1:solver.ny-1;
        v_ly = -solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2 * ones(1,solver.ny-1);

        solver.substrate_pars(si).My = sparse([r_ind_uy,r_ind_dy,r_ind_ly],...
            [c_ind_uy,c_ind_dy,c_ind_ly],...
            [v_uy    ,v_dy   ,v_ly    ]);

        %% set up in z direction
        r_ind_uz = 1:solver.nz-1;
        c_ind_uz = 2:solver.nz;
        v_uz = -solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2 * ones(1,solver.nz-1);

        r_ind_dz = 1:solver.nz;
        c_ind_dz = 1:solver.nz;
        v_dz(1,2:solver.nz-1) = 1 + solver.dt*solver.substrate_pars(si).degradation/3 + 2*solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2;
        v_dz([1,solver.nz]) = 1 + solver.dt*solver.substrate_pars(si).degradation/3 + solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2;

        r_ind_lz = 2:solver.nz;
        c_ind_lz = 1:solver.nz-1;
        v_lz = -solver.dt*solver.substrate_pars(si).diffusion/solver.cell_width^2 * ones(1,solver.nz-1);

        solver.substrate_pars(si).Mz = sparse([r_ind_uz,r_ind_dz,r_ind_lz],...
            [c_ind_uz,c_ind_dz,c_ind_lz],...
            [v_uz    ,v_dz   ,v_lz    ]);

    end

elseif solver.method=="global"

    solver.dt = solver.Tf;
    solver.M = accumarray([reshape(solver.regions(1:end-1,:,:),[],1),reshape(solver.regions(2:end,:,:),[],1);... % count neighbors to left
                    reshape(solver.regions(2:end,:,:),[],1),reshape(solver.regions(1:end-1,:,:),[],1);... % count neighbors to right
                    reshape(solver.regions(:,1:end-1,:),[],1),reshape(solver.regions(:,2:end,:),[],1);... % count neighbors in front
                    reshape(solver.regions(:,2:end,:),[],1),reshape(solver.regions(:,1:end-1,:),[],1);... % count neighbors behind
                    reshape(solver.regions(:,:,1:end-1),[],1),reshape(solver.regions(:,:,2:end),[],1);... % count neighbors below
                    reshape(solver.regions(:,:,2:end),[],1),reshape(solver.regions(:,:,1:end-1),[],1);... % count neighbors above 
                    reshape(solver.regions(1,:,:),[],1),reshape(solver.regions(1,:,:),[],1);... % count neighbors left of left face as same region as center
                    reshape(solver.regions(end,:,:),[],1),reshape(solver.regions(end,:,:),[],1);... % count neighbors right of right face as same region as center
                    reshape(solver.regions(:,1,:),[],1),reshape(solver.regions(:,1,:),[],1);... % count neighbors in front of front face as same region as center
                    reshape(solver.regions(:,end,:),[],1),reshape(solver.regions(:,end,:),[],1);... % count neighbors behind back face as same region as center
                    reshape(solver.regions(:,:,1),[],1),reshape(solver.regions(:,:,1),[],1);... % count neighbors below bottom face as same region as center
                    reshape(solver.regions(:,:,end),[],1),reshape(solver.regions(:,:,end),[],1);... % count neighbors above top face as same region as center
                    ],1);
    assert(issymmetric(solver.M))
    solver.M = solver.M./sum(solver.M,1);
    solver.M = solver.M/solver.cell_width^2;

    solver.region_volumes = accumarray(solver.regions(:),1);
    solver.regional_BV_prop = (accumarray(solver.regions(:),solver.IE(:))./solver.region_volumes)';

    solver.num_cell_bound = length(solver.bound_inds)*solver.num_types*solver.n_regions;
    [~,solver.perm_order] = sort([solver.bound_inds,solver.substrate_inds]);


else
    error('unknown method')
end

for si = 1:solver.num_substrates
    %% entry locations of inhibitor
    if solver.is_pk(si)

        % x' = f*(y-x); % rate of change of concentration in perivascular region of TME
        % y' = k12*(z-y) - d*y; % rate of change of concentration in circulation
        % z' = k21*(y-z); % rate of change of concentration in periphery (as separate compartment from TME)

        pk_matrix = zeros(3); % row 1: concentration in microenvironment at blood vessel; row 2: concentration in blood vessel; row 3: concentration in periphery
        pk_matrix(1,1:2) = solver.substrate_pars(si).fluid_exchange_rate*[-1,1]; % update concentration in TME by blood vessel based on difference with circulation concentration
        pk_matrix(2,2:3) = solver.substrate_pars(si).k12*[-1,1]; % intercompartmental clearance
        pk_matrix(2,2) = pk_matrix(2,2) - solver.substrate_pars(si).sysdecay; % systemic clearance
        pk_matrix(3,2:3) = solver.substrate_pars(si).k21*[1,-1]; % intercompartmental clearance

        [solver.substrate_pars(si).eV,solver.substrate_pars(si).ev] = eig(pk_matrix);
        solver.substrate_pars(si).eVinv = inv(solver.substrate_pars(si).eV);

        if solver.method=="local"
            solver.substrate_pars(si).perivascular_update_dual =  solver.substrate_pars(si).eV(1,:)*(exp(solver.dt*diag(solver.substrate_pars(si).ev)).*solver.substrate_pars(si).eVinv);
            solver.substrate_pars(si).circulation_update_dual =  solver.substrate_pars(si).eV(2,:)*(exp(solver.dt*diag(solver.substrate_pars(si).ev)).*solver.substrate_pars(si).eVinv);
        elseif solver.method=="global"
            solver.substrate_pars(si).circulation_update_dual = @(dt) solver.substrate_pars(si).eV(2,:)*(exp(dt*diag(solver.substrate_pars(si).ev)).*solver.substrate_pars(si).eVinv);
        end
        solver.substrate_pars(si).periphery_update_dual = solver.substrate_pars(si).eV(3,:)*(exp(solver.dt*diag(solver.substrate_pars(si).ev)).*solver.substrate_pars(si).eVinv);
    end
end

solver.setup_done = true;
