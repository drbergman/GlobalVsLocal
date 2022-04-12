clearvars;

% This script will test the effects of changing PK parameters on the
% equivalence of our two methods

base_name = 'varyPK';
nsamps = 1;
min_parfor_num = 4; % since these runs are not being timed, can possibly use a parpool to speed up
N0 = 1e4; % initial number of tumor cells
save_output = false; % set to true if you want to save the output

%% dosing regimens
aFGFR3_start_min = 1; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1; % days between successive doses of anti-FGFR3
DoW_start = 6; % day of week to start on; 0 = monday (so 6 means the first simmed day is considered a Sunday and the first dose should be given on Day 1 (see aFGFR3_start_simday below)
censor_date = 8; % day in simulation to end

n_doses_aFGFR3 = 5; % max number of doses of FGFR3 to give

%% initialize pars
pars = basePars_FGFR3();
pars.main.max_pde_dt = (.1/60)/24; % smaller pde time steps needed for faster diffusion and degradation
aFGFR3_diffusion = 1e-5*logspace(-2,0,3);
aFGFR3_circ0 = 1e3*logspace(-1,1,3);
aFGFR3_degradation = 144*logspace(-2,0,3);

aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min; % simulation day to start therapy
events = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]);

f = @(pars) startPatient(N0,pars,events{1}); % function to run simulation with a given set of parameters

method = {'Local','Global'};
sz = [length(method),nsamps,length(aFGFR3_diffusion),length(aFGFR3_circ0),length(aFGFR3_degradation)];

times = zeros(sz);
total_runs = prod(sz);

if total_runs>=min_parfor_num
    F(1:total_runs) = parallel.FevalFuture;
end

mu_n = 0; % weighted average of simulation times for displaying estimated time remaining for cohort
timer_start = tic;

for i = total_runs:-1:1
    [mi,si,adiffi,acirci,adegi] = ind2sub(sz,i);
    q = pars;
    q.main.method = method{mi};
    q.main.aFGFR3_diffusion = aFGFR3_diffusion(adiffi);
    q.main.aFGFR3_circ0 = aFGFR3_circ0(acirci);
    q.main.aFGFR3_degradation = aFGFR3_degradation(adegi);
    
    if total_runs<min_parfor_num
        timerVal = tic;
        TRACKED(mi,si,adiffi,acirci,adegi) = f(q);
        times(mi,si,adiffi,acirci,adegi) = toc(timerVal);
        
        % displays estimated time remaining for whole cohort (not part of the
        % method, but helpful to let you know how long to expect)
        num_done = total_runs - i + 1;
        mu_n = ((num_done-1)*mu_n+2*times(mi,si,adiffi,acirci,adegi))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
        etr = mu_n * (total_runs-num_done);
        fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
            num_done,total_runs,100*num_done/total_runs,duration(0,0,times(mi,si,adiffi,acirci,adegi)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
       
    else
        F(i) = parfeval(f,1,q);
    end
end

if total_runs >= min_parfor_num
    this_parpool = gcp;
    numworkers = this_parpool.NumWorkers;
    for i = 1:total_runs
        if mod(i,4)==1
            tic;
        end
        [idx,next_T] = fetchNext(F);
        [mi,si,adiffi,acirci,adegi] = ind2sub(sz,idx);
        TRACKED(mi,si,adiffi,acirci,adegi) = next_T;
        if mod(i,numworkers)==0 % update the estimated time remaining after all workers have returned
            v_n = toc;
            mu_n = (((i/4)-1)*mu_n+2*v_n)/((i/4)+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
            etr = mu_n*(total_runs-i)/4;
            fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
                i,total_runs,100*i/total_runs,duration(0,0,v_n),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        end
    end
    
end

if save_output
    if ~exist('./data','dir') % make sure the data directory exists
        mkdir('./data')
    end
    filename = nextFileName('data',base_name,3);
    save(filename,'-v7.3')
end
