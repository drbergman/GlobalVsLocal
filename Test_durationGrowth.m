close all; clear all; 
% This script will test how the speed of the models grows with N0

base_name = 'durationGrowth';

n = 7;
nsamps = 1;

%% parameters that are constant for all patients
N0 = round(logspace(log10(50),log10(50000),n)); % initial number of tumor cells

%% dosing regimens
aFGFR3_start_min = 1; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1;
DoW_start = 6; % day of week to start on; 0 = monday
censor_date = 8;

n_doses_aFGFR3 = 10;


%% set file names
filename = nextFileName('data',base_name,3);

%% initialize pars
pars = basePars_FGFR3();
pars.main.max_pde_dt = (.25/60)/24;


aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min;
EVENTS = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]);
n_subcohorts = numel(EVENTS);

f = @(ni,pars) startPatient(N0(ni),pars,EVENTS{1});

method = {'Local','Global'};

sz = [length(method),nsamps,n];
TIMES = zeros(sz);
total_runs = prod(sz);

mu_n = 0;
timer_start = tic;
            
for i = total_runs:-1:1
    [mi,si,ni] = ind2sub(sz,i);
    q = pars;
    q.main.method = method{mi};
    
    timerVal = tic;
    TRACKED(mi,si,ni) = f(ni,q);
    TIMES(mi,si,ni) = toc(timerVal);
    
    num_done = total_runs - i + 1;
    mu_n = ((num_done-1)*mu_n+2*TIMES(mi,si,ni))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
    etr = mu_n * (total_runs-num_done);
    fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
        num_done,total_runs,100*num_done/total_runs,duration(0,0,TIMES(mi,si,ni)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        
end

%% clear variables I definitely do not want to save

% save(filename,'-v7.3')
