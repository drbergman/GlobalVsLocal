close all; clear all;
% This script will test the effects of changing PK parameters on the
% equivalence of our two methods

base_name = 'varyReactionRates';

nsamps = 1;
min_parfor_num = 4;

%% parameters that are constant for all patients
N0 = 1e4; % initial number of tumor cells

%% dosing regimens
aFGFR3_start_min = 1; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1;
DoW_start = 6; % day of week to start on; 0 = monday
censor_date = 8;

n_doses_aFGFR3 = 100;


%% set file names
filename = nextFileName('data',base_name,3);

%% initialize pars
pars = basePars_FGFR3();
pars.main.max_pde_dt = (.1/60)/24;

k_on_R = 1e2*logspace(-1,1,3); % used for both monomers and dimers
k_off_R = 7.8e2*logspace(-1,1,3); % used for both monomers and dimers
kf = 460.8*logspace(-1,1,3);
kr = 400*logspace(-1,1,3);

aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min;
EVENTS = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]);
n_subcohorts = numel(EVENTS);

f = @(pars) startPatient(N0,pars,EVENTS{1});

method = {'Local','Global'};

sz = [length(method),nsamps,length(k_on_R),length(k_off_R),length(kf),length(kr)];
method_ind = 1;
sample_ind = 2;
k_on_R_ind = 3;
k_off_R_ind = 4;
kf_ind = 5;
kr_ind = 6;

TIMES = zeros(sz);
total_runs = prod(sz);

if total_runs>=min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
end

mu_n = 0;
timer_start = tic;

for i = total_runs:-1:1
    [method_ind,si,k_on_Ri,k_off_Ri,kfi,kri] = ind2sub(sz,i);
    q = pars;
    q.main.method = method{method_ind};
    q.main.k_on_R = k_on_R(k_on_Ri);
    q.main.k_on_D = k_on_R(k_on_Ri);
    q.main.k_off_R = k_off_R(k_off_Ri);
    q.main.k_off_D = k_off_R(k_off_Ri);
    
    q.main.kf = kf(kfi);
    q.main.kr = kr(kri);
    
    if total_runs<min_parfor_num
        timerVal = tic;
        TRACKED(method_ind,si,k_on_Ri,k_off_Ri,kfi,kri) = f(q);
        TIMES(method_ind,si,k_on_Ri,k_off_Ri,kfi,kri) = toc(timerVal);
        
        num_done = total_runs - i + 1;
        mu_n = ((num_done-1)*mu_n+2*TIMES(method_ind,si,k_on_Ri,k_off_Ri,kfi,kri))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
        etr = mu_n * (total_runs-num_done);
        fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
            num_done,total_runs,100*num_done/total_runs,duration(0,0,TIMES(method_ind,si,k_on_Ri,k_off_Ri,kfi,kri)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        
    else
        F(i) = parfeval(f,1,q);
    end
    
end

if total_runs >= min_parfor_num
    for i = 1:total_runs
        if mod(i,4)==1
            tic;
        end
        [idx,next_T] = fetchNext(F);
        [method_ind,si,k_on_Ri,k_off_Ri,kfi,kri] = ind2sub(sz,idx);
        TRACKED(method_ind,si,k_on_Ri,k_off_Ri,kfi,kri) = next_T;
        if mod(i,4)==0
            v_n = toc;
            mu_n = (((i/4)-1)*mu_n+2*v_n)/((i/4)+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
            etr = mu_n*(total_runs-i)/4;
            fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
                i,total_runs,100*i/total_runs,duration(0,0,v_n),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        end
    end
    
end

%% clear variables I definitely do not want to save
% clear F next_T
% 
% save(filename,'-v7.3')
