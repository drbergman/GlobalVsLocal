close all; clear all;
% This script will test what happens when the parameters are changed so
% that the inhibitor behaves like an antibody

base_name = 'AbDynamics';

nsamps = 1;
min_parfor_num = 4;

%% parameters that are constant for all patients
N0 = 1e4; % initial number of tumor cells

%% dosing regimens
aFGFR3_start_min = 1; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 3;
DoW_start = 6; % day of week to start on; 0 = monday
censor_date = 20;

n_doses_aFGFR3 = 100;


%% set file names
filename = nextFileName('data',base_name,3);

%% initialize pars
pars = basePars_FGFR3();
pars.main.max_pde_dt = (.1/60)/24;

pars.main.aFGFR3_influx = 1;
pars.main.aFGFR3_diffusion = 1e-7;
pars.main.aFGFR3_circ0 = 1e0;
pars.main.aFGFR3_degradation = 1.44;
pars.main.aFGFR3_sysdecay = log(2)/1;

pars.main.k_on_R = 100e4 * 1e-9 * 3600 * 24; % used for both monomers and dimers
pars.main.k_off_R = 1e-4 * 3600 * 24; % used for both monomers and dimers
pars.main.k_on_D = 100e4 * 1e-9 * 3600 * 24; % used for both monomers and dimers
pars.main.k_off_D = 1e-4 * 3600 * 24; % used for both monomers and dimers
pars.main.kf = 460.8;
pars.main.kr = 400;

influx_method = {'everywhere','line','outside'};


aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min;
EVENTS = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]);
n_subcohorts = numel(EVENTS);

f = @(pars) startPatient(N0,pars,EVENTS{1});

method = {'Local','Global'};
figure(2);
for i = 1:length(influx_method)
    subplot(length(influx_method),1,i); hold on;
    title(influx_method{i})
end

color = parula(2);

sz = [length(method),nsamps,length(influx_method)];
sample_ind = 2;
method_ind = 1;
influx_ind = 3;

TIMES = zeros(sz);
total_runs = prod(sz);

if total_runs>=min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
end

mu_n = 0;
timer_start = tic;

for i = total_runs:-1:1
    [method_ind,si,influx_method_ind] = ind2sub(sz,i);
    q = pars;
    q.main.method = method{method_ind};
    q.main.inhibitor_entry = influx_method{influx_method_ind};
    
    switch influx_method{influx_method_ind} % empirically, there's less inhibitor in other methods, this will try to get them all to have the same amount of inhibitor in the TME
        case 'line'
            q.main.aFGFR3_circ0 = q.main.aFGFR3_circ0*460.544/.0656948;
        case 'outside'
            q.main.aFGFR3_circ0 = q.main.aFGFR3_circ0*460.544/372.886;
    end
    
    if total_runs<min_parfor_num
        timerVal = tic;
        TRACKED(method_ind,si,influx_method_ind) = f(q);
        TIMES(method_ind,si,influx_method_ind) = toc(timerVal);
        
        num_done = total_runs - i + 1;
        mu_n = ((num_done-1)*mu_n+2*TIMES(method_ind,si,influx_method_ind))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
        etr = mu_n * (total_runs-num_done);
        fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
            num_done,total_runs,100*num_done/total_runs,duration(0,0,TIMES(method_ind,si,influx_method_ind)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        
        subplot(length(influx_method),1,influx_method_ind)
        plot(TRACKED(method_ind,si,influx_method_ind).T,TRACKED(method_ind,si,influx_method_ind).NT,'Color',color(method_ind,:))
        drawnow
        
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
        [method_ind,si,influx_method_ind] = ind2sub(sz,idx);
        TRACKED(method_ind,si,influx_method_ind) = next_T;
        if mod(i,4)==0
            v_n = toc;
            mu_n = (((i/4)-1)*mu_n+2*v_n)/((i/4)+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
            etr = mu_n*(total_runs-i)/4;
            fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
                i,total_runs,100*i/total_runs,duration(0,0,v_n),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        end
        subplot(length(influx_method),1,influx_method_ind)
        plot(TRACKED(method_ind,si,influx_method_ind).T,TRACKED(method_ind,si,influx_method_ind).NT,'Color',color(method_ind,:))
        drawnow
    end
    
end

%% clear variables I definitely do not want to save
% clear F next_T
% 
% save(filename,'-v7.3')
