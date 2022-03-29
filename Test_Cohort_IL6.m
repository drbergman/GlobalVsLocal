clearvars;
close all
%% set up cohort
nsamps = 1;
method = ["local","global"];
use_pde_in_global = false; % if not toggling pde in global method, which version of global method to use
base_name = "il6_cohort";
min_parfor = 4;
%% set up common parameters
input = basePars_IL6();

%% some parameters you may want to change for this cohort

% input.flags.plotFigs = true;
% input.plot_properties.plotLocations = true;
% input.solver.max_pde_dt = 1;
% input.pars.desired_dt = 1;
% input.update_pars.P_ec50 = 2;
% input.pars.TME_size = [20,20,20];
% input.simpars.censor_date = 0.1;
% input.solver.dirichlet_condition(1) = 5e3;
% input.initialization_pars.N0 = [100,10,20];
% input.solver.substrate_pars(1).degradation = 0.4;
% input.solver.substrate_pars(1).diffusion = 750;
% input.solver.plot_every_pde_step = true;
% input.event_pars.prolif_rate = input.event_pars.prolif_rate;
% input.event_pars.delta = 5*input.event_pars.delta;
% input.update_pars.mu_S = 0.6;
% input.event_pars.move_rate_in_microns = 0.2;

% input.solver.agent_ode_pars.kf = 0; % dont let cells affect il6 concentrations

% input.initialization_pars.start_with_drug(2) = true; % start with aIL6R
% input.initialization_pars.circulation_concentration(2) = 1e2;
% input.solver.substrate_pars(2).sysdecay = 0.04/(60*24); % set so half-life in circulation is log(2)/0.04 ~= 17 days, which we can see in one simulation
% input.solver.substrate_pars(2).fluid_exchange_rate = 20*input.solver.substrate_pars(2).k12; % since 1/20 of the TME (in a 20x20x20 TME) exchanges with circulation, this rate should mean this part of the periphery gets as much aIL6R as the rest of the periphery

%% run cohort

sz = [length(method),nsamps];
total_runs = prod(sz);
F(1:total_runs) = parallel.FevalFuture;

for mi = length(method):-1:1
    input.method = method(mi);
    for si = nsamps:-1:1
        if total_runs>=min_parfor
            ind = sub2ind(sz,mi,si);
            F(ind) = parfeval(@simPatient_IL6,1,input);
        else
            out(si,mi) = simPatient_IL6(input);
        end
    end
end

if total_runs>=min_parfor
    for i = 1:total_runs
        [idx,next_out] = fetchNext(F);
        [mi,si] = ind2sub(sz,idx);
        out(si,mi) = next_out;
    end
end

%% save output
% filename = nextFileName('data',base_name,3);
% clear F next_out
% save(filename,'-v7.3')
