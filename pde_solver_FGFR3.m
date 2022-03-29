function [aFGFR3,aFGFR3_circ,p,tumors,timer] = pde_solver_FGFR3(xx,yy,zz,aFGFR3,...
    aFGFR3_circ,p,tumors,tumor_receptors_inds,ind_ind,timer,kf,kr,kp,...
    k_on_R,k_off_R,k_on_D,k_off_D,inhibitor_entry,location_inds,cell_width,vessel_spacing)

nx = length(xx);
ny = length(yy);
nz = length(zz);

dx = xx(2)-xx(1);
dy = yy(2)-yy(1);
dz = zz(2)-zz(1);

tt = linspace(0,p.t_pde,p.nt);
dt = tt(2)-tt(1);

IC = tumors(:,tumor_receptors_inds);
NT = size(tumors,1);

if ~p.setup_done % then need to setup solver stuff
    %% set up in x direction
    r_ind_ux = 1:nx-1;
    c_ind_ux = 2:nx;
    v_ux = -dt*p.aFGFR3_diffusion/dx^2 * ones(1,nx-1);
    
    r_ind_dx = 1:nx;
    c_ind_dx = 1:nx;
    v_dx(1,2:nx-1) = 1 + dt*p.aFGFR3_degradation/3 + 2*dt*p.aFGFR3_diffusion/dx^2;
    v_dx([1,nx]) = 1 + dt*p.aFGFR3_degradation/3 + dt*p.aFGFR3_diffusion/dx^2;
    
    r_ind_lx = 2:nx;
    c_ind_lx = 1:nx-1;
    v_lx = -dt*p.aFGFR3_diffusion/dx^2 * ones(1,nx-1);
    
    p.Mx = sparse([r_ind_ux,r_ind_dx,r_ind_lx],...
        [c_ind_ux,c_ind_dx,c_ind_lx],...
        [v_ux    ,v_dx   ,v_lx     ],nx,nx);
    
    %% set up in y direction
    r_ind_uy = 1:ny-1;
    c_ind_uy = 2:ny;
    v_uy = -dt*p.aFGFR3_diffusion/dy^2 * ones(1,ny-1);
    
    r_ind_dy = 1:ny;
    c_ind_dy = 1:ny;
    v_dy(1,2:ny-1) = 1 + dt*p.aFGFR3_degradation/3 + 2*dt*p.aFGFR3_diffusion/dy^2;
    v_dy([1,ny]) = 1 + dt*p.aFGFR3_degradation/3 + dt*p.aFGFR3_diffusion/dy^2;
    
    r_ind_ly = 2:ny;
    c_ind_ly = 1:ny-1;
    v_ly = -dt*p.aFGFR3_diffusion/dy^2 * ones(1,ny-1);
    
    p.My = sparse([r_ind_uy,r_ind_dy,r_ind_ly],...
        [c_ind_uy,c_ind_dy,c_ind_ly],...
        [v_uy    ,v_dy   ,v_ly    ]);
    
    %% set up in z direction
    r_ind_uz = 1:nz-1;
    c_ind_uz = 2:nz;
    v_uz = -dt*p.aFGFR3_diffusion/dz^2 * ones(1,nz-1);
    
    r_ind_dz = 1:nz;
    c_ind_dz = 1:nz;
    v_dz(1,2:nz-1) = 1 + dt*p.aFGFR3_degradation/3 + 2*dt*p.aFGFR3_diffusion/dz^2;
    v_dz([1,nz]) = 1 + dt*p.aFGFR3_degradation/3 + dt*p.aFGFR3_diffusion/dz^2;
    
    r_ind_lz = 2:nz;
    c_ind_lz = 1:nz-1;
    v_lz = -dt*p.aFGFR3_diffusion/dz^2 * ones(1,nz-1);
    
    p.Mz = sparse([r_ind_uz,r_ind_dz,r_ind_lz],...
        [c_ind_uz,c_ind_dz,c_ind_lz],...
        [v_uz    ,v_dz   ,v_lz    ]);
    
    %% entry locations of inhibitor
    switch inhibitor_entry
        case 'everywhere'
            p.IE = true(size(aFGFR3));
        case 'border'
            p.IE = true(size(aFGFR3));
            p.IE(2:end-1,2:end-1,2:end-1) = false;
        case 'line'
            p.IE = false(size(aFGFR3));
            [~,xind] = min(abs(xx));
            [~,yind] = min(abs(yy));
            p.IE(xind,yind,:) = true;
        case 'plane'
            p.IE = false(size(aFGFR3));
            [~,xind] = min(abs(xx));
            p.IE(xind,:,:) = true;
        case 'outside'
            % this is handled below
        case 'multiple lines'
            p.IE = false(size(aFGFR3));
            xloc_bv = vessel_spacing*(ceil(xx(1)/vessel_spacing):floor(xx(end)/vessel_spacing));
            yloc_bv = vessel_spacing*(ceil(yy(1)/vessel_spacing):floor(yy(end)/vessel_spacing));
            [~,xind_bv] = min(abs(xx-xloc_bv'),[],2);
            [~,yind_bv] = min(abs(yy-yloc_bv'),[],2);
            p.IE(xind_bv,yind_bv,:) = true;
        otherwise
            error('Have not added ''%s'' to the local method\n',inhibitor_entry)
    end
    
    p.setup_done = true;
end

Js = zeros(6,6,NT);
Js(1,2,:) = 2*(kr+kp);
Js(1,4,:) = k_off_R;
Js(2,6,:) = k_off_D;
Js(3,4,:) = k_off_R;
Js(3,6,:) = k_off_D;
Js(4,5,:) = 2*(kr+kp);
Js(5,5,:) = -kr-kp;
Js(6,6,:) = -k_off_D;

if strcmp(inhibitor_entry,'outside')
    com = cell_width*mean(tumors(:,location_inds),1); % tumor center of mass
    max_distance = cell_width^2*max(sum((tumors(:,location_inds)-com).^2,2)); % technically it's the square of the distance
    p.IE = find(sum((allCombos(xx,yy,zz,'matlab')-com).^2,2)>max_distance);
end

for ti = 2:p.nt
    
    timer.pde = tic;
    
    %% diffusion decay terms
    % x-diffusion
    for j = 1:ny
        for k = 1:nz
            aFGFR3(:,j,k) = p.Mx\aFGFR3(:,j,k);
        end
    end
    
    % y-diffusion
    for i = 1:nx
        for k = 1:nz
            aFGFR3(i,:,k) = p.My\reshape(aFGFR3(i,:,k),[],1);
        end
    end
    
    % z-diffusion
    for i = 1:nx
        for j = 1:ny
            aFGFR3(i,j,:) = p.Mz\reshape(aFGFR3(i,j,:),[],1);
        end
    end
    
    %% Bulk supply/uptake
    aFGFR3(p.IE) = aFGFR3(p.IE)*exp(-p.aFGFR3_eflux*dt) + p.aFGFR3_influx*aFGFR3_circ*(exp(-p.aFGFR3_sysdecay*dt)-exp(-(p.aFGFR3_eflux+p.aFGFR3_degradation)*dt))/(p.aFGFR3_eflux+p.aFGFR3_degradation - p.aFGFR3_sysdecay); % this is the exact value when just considering influx and eflux between periphery and circulation and also counting the degradation of the newly-influxed inhibitor
    
    timer.pde_duration = timer.pde_duration + toc(timer.pde);
    
    %% Cell supply/uptake
    IC(:,3) = aFGFR3(tumors(:,ind_ind));
    timer.ode = tic;
    dIC = agentODE_FGFR3(IC,kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D); % working on this here now
    Js = agentODEJacobian_FGFR3(Js,IC,kf,kr,kp,k_on_R,k_off_R,k_on_D); % .5 to say (x2-x1)/dt = (f(x1)+f(x2))/2 OR (x2-x1)/dt = f((x1+x2)/2); this .5 is now accounted for in constructing Js
    for ci = 1:NT
        IC(ci,:) = firstOrderImplicitUpdate(dt,IC(ci,:)',dIC(ci,:)',Js(:,:,ci),kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D);
    end
    
    aFGFR3(tumors(:,ind_ind)) = IC(:,3); % update locations that have a tumor cell
    
    timer.ode_duration = timer.ode_duration+toc(timer.ode);
    
    %% update circulating inhibitor
    
    aFGFR3_circ = aFGFR3_circ*exp(-p.aFGFR3_sysdecay*dt);
    
end
tumors(:,tumor_receptors_inds) = IC;

end

function x = firstOrderImplicitUpdate(dt,IC,dIC,J,kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D)
x = (eye(6)-.5*dt*J)\(IC + dt*(dIC-.5*J*IC));
if any(x<0)
    x = IC;
    for i = 1:2
        x = firstOrderImplicitUpdate(dt/2,x,dIC,J,kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D);
        if i==1
            dIC = agentODE_FGFR3(x',kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D)';
            J = agentODEJacobian_FGFR3(J,x',kf,kr,kp,k_on_R,k_off_R,k_on_D);
        end
    end
end
end