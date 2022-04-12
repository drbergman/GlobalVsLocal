clearvars;

% This script will test the effects of changing the location of influx of 
% drug on the equivalence of our two methods

base_name = 'varyInflux';
nsamps = 1;
min_parfor_num = 4;
N0 = 1e4; % initial number of tumor cells
save_output = false; % set to true if you want to save the output

%% dosing regimens
aFGFR3_start_min = 1; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1;
DoW_start = 6; % day of week to start on; 0 = monday
censor_date = 8;

n_doses_aFGFR3 = 100;

%% initialize pars
pars = basePars_FGFR3();
pars.main.max_pde_dt = (.1/60)/24;
influx_method = {'everywhere','line','outside'};

aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min;
events = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]);

f = @(pars) startPatient(N0,pars,events{1});

method = {'Local','Global'};

sz = [length(method),nsamps,length(influx_method)];

times = zeros(sz);
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
            q.main.aFGFR3_circ0 = q.main.aFGFR3_circ0*61.1866/0.0114637;
        case 'outside'
            q.main.aFGFR3_circ0 = q.main.aFGFR3_circ0*61.1866/53.2994;
    end
    if total_runs<min_parfor_num
        timerVal = tic;
        TRACKED(method_ind,si,influx_method_ind) = f(q);
        times(method_ind,si,influx_method_ind) = toc(timerVal);
        
        num_done = total_runs - i + 1;
        mu_n = ((num_done-1)*mu_n+2*times(method_ind,si,influx_method_ind))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
        etr = mu_n * (total_runs-num_done);
        fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
            num_done,total_runs,100*num_done/total_runs,duration(0,0,times(method_ind,si,influx_method_ind)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
        
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
    end
    
end

if save_output
    if ~exist('./data','dir') % make sure the data directory exists
        mkdir('./data')
    end
    filename = nextFileName('data',base_name,3);
    save(filename,'-v7.3')
end
