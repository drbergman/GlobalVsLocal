close all; clear all; 
% This script will test the speed of our model running local vs global and
% their equivalence

base_name = 'ourModel';

nsamps = 1;

%% parameters that are constant for all patients
N0 = 1e4; % initial number of tumor cells

%% dosing regimens
aFGFR3_start_min = 0; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1;
DoW_start = 6; % day of week to start on; 0 = monday
censor_date = 8;

n_doses_aFGFR3 = 15;

%% set file names
filename = nextFileName('data',base_name,3);

%% initialize pars
pars = basePars_FGFR3();
pars.main.dt = 9/24;
pars.main.max_pde_dt = (10/60)/24;
aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min;
EVENTS = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]);
n_subcohorts = numel(EVENTS);

f = @(pars) startPatient_FGFR3(N0,pars,EVENTS{1});

method = {'Local','Global'};

method_ind = 1;
sample_ind = 2;
sz = [length(method),nsamps];
total_runs = prod(sz);
TIMES = zeros(sz);

mu_n = 0; % weighted average of simulation times for displaying estimated time remaining for cohort
timer_start = tic;

for i = total_runs:-1:1
    [mi,si] = ind2sub(sz,i);
    q = pars;
    q.main.method = method{mi};
    
    timerVal = tic;
    TRACKED(mi,si) = f(q);
    TIMES(mi,si) = toc(timerVal);
    
    % displays estimated time remaining for whole cohort
    num_done = total_runs - i + 1;
    mu_n = ((num_done-1)*mu_n+2*TIMES(mi,si))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
    etr = mu_n * (total_runs-num_done);
    fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
        num_done,total_runs,100*num_done/total_runs,duration(0,0,TIMES(mi,si)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
end

%% clear variables I definitely do not want to save

% save(filename,'-v7.3')
