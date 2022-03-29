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

% input.pars.desired_dt = 1;
% input.pars.TME_size = [20,20,20];
% input.simpars.censor_date = 0.1;
% input.initialization_pars.N0 = [100,10,20];

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