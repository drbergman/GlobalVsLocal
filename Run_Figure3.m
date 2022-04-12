clearvars;

% This script will test the speed of our model running local vs global and
% their equivalence. it generates the data necessary for Figure 3.

base_name = 'Figure2_';
nsamps = 3; % number of samples to run for each method
N0 = 1e1; % initial number of tumor cells
save_output = false; % set to true if you want to save the output

%% dosing regimens
aFGFR3_start_min = 0; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1; % days between successive doses of anti-FGFR3
DoW_start = 6; % day of week to start on; 0 = monday (so 6 means the first simmed day is considered a Sunday and the first dose should be given on Day 1 (see aFGFR3_start_simday below)
censor_date = 2; % day in simulation to end

n_doses_aFGFR3 = 5; % max number of doses of FGFR3 to give

%% set filename for saving (if desired)

%% initialize pars
pars = basePars_FGFR3();

aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min; % simulation day to start therapy

events = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]); % decide on when events occur within each simulation; events include new doses and censoring

f = @(pars) startPatient_FGFR3(N0,pars,events{1}); % function to run simulation with a given set of parameters

method = {'Local','Global'};

% organizing the output of the cohort
method_ind = 1;
sample_ind = 2;
sz = [length(method),nsamps];

% preparing the timing of the runs
total_runs = prod(sz);
times = zeros(sz);

mu_n = 0; % weighted average of simulation times for displaying estimated time remaining for cohort
timer_start = tic;

for i = total_runs:-1:1
    [mi,si] = ind2sub(sz,i); % first loop over methods to get better approximation of time remaining earlier on (see below)
    q = pars;
    q.main.method = method{mi}; % set the method to be used
    
    timerVal = tic;
    tracked(mi,si) = f(q); % run (and time) the simulation
    times(mi,si) = toc(timerVal);
    
    % displays estimated time remaining for whole cohort (not part of the
    % method, but helpful to let you know how long to expect)
    num_done = total_runs - i + 1;
    mu_n = ((num_done-1)*mu_n+2*times(mi,si))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
    etr = mu_n * (total_runs-num_done);
    fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
        num_done,total_runs,100*num_done/total_runs,duration(0,0,times(mi,si)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))
end

%% clear variables I definitely do not want to save

if save_output
    if ~exist('./data','dir') % make sure the data directory exists
        mkdir('./data')
    end
    filename = nextFileName('data',base_name,3);
    save(filename,'-v7.3')
end

%% generate panels

nfigs = 6;
figs = gobjects(nfigs,1);
method_colors = lines(2);

t = tracked(1).T;

%% tumor cell counts plots
all_tumor_sizes = zeros(length(t),nsamps,2);
for mi = 1:2
    for si = 1:nsamps
        all_tumor_sizes(:,si,mi) = tracked(mi,si).NT;
    end
end

figs(1) = figure("Name","Tumor Growth Curves");
hold on
for mi = 1:2
    plot(t,all_tumor_sizes(:,:,mi),'Color',method_colors(mi,:))
end
xlabel('Time (d)')
ylabel('Tumor Cell Count')

figs(2) = figure("Name","Residual Tumor Growth Curves");
hold on;
for mi = 1:2
    plot(t,all_tumor_sizes(:,:,mi) - mean(all_tumor_sizes,2:3),"Color",method_colors(mi,:))
end
xlabel('Time (d)')
ylabel('Residual Tumor Cell Count')

%% expected growth plots
all_eg = zeros(length(t)-1,nsamps,2); % note this is not computed at the final time point
for mi = 1:2
    for si = 1:nsamps
        all_eg(:,si,mi) = tracked(mi,si).simple_expectedGrowth; % the field expected_growth is the stochastic expected growth shown in the SI
    end
end

figs(3) = figure("Name","Expected Growth Rate");
hold on
for mi = 1:2
    plot(t(1:end-1),all_eg(:,:,mi),"Color",method_colors(mi,:))
end
xlabel('Time (d)')
ylabel('Expected Growth Rate (d^{-1})')

% note that this computation is made by simulation rather than by mean.
% This shows there is little variation between simulations in the average
% signaling in a tumor
figs(4) = figure("Name","Difference in Expected Growth Rate");
plot(t(1:end-1),diff(all_eg,1,3),"Color",method_colors(mi,:))
xlabel('Time (d)')
ylabel('Expected Growth Rate (d^{-1})')

%% percent difference in tumor size
int_eg = pars.main.dt*cumsum(cat(1,zeros(1,nsamps),diff(all_eg,1,3))); % integrate all the expected growth rates and compute the difference between the local and global methods
percent_diff_tegs = mean(100*(exp(int_eg)-1),2); % exponentiate (difference of) TEG to get (ratio of) expected tumor size; then convert to percent difference x-->100*(x-1) and compute mean

all_tumor_sizes_diffs = diff(all_tumor_sizes,1,3);
percent_diffs = 100*all_tumor_sizes_diffs./all_tumor_sizes(:,:,1);

figs(5) = figure("Name","Percent Difference in Total Expected Growth");
hold on
plot(t,percent_diffs,'Color',.7*ones(1,3));
plot(t,percent_diff_tegs,'Color',[.2,.1,.8],'LineWidth',3,'DisplayName','Difference Using TEG');
xlabel('Time (d)')
ylabel('Percent Difference')

%% total wall times
figs(6) = figure("Name","Total Simulation Wall Times");
hold on
for mi = 1:2
    histogram(log10(times(mi,:)),'FaceColor',method_colors(mi,:),'Normalization','count','EdgeColor','none','FaceAlpha',1);
end
xL = [floor(log10(min(times,[],'all'))),ceil(log10(max(times,[],'all')))];
xlim(xL)
xt = ceil(xL(1)):floor(xL(2));
xticks(xt)
xticklabels(10.^xt)
xlabel('Wall Time (s)')
ylabel('Frequency')
legend('Local','Global')




