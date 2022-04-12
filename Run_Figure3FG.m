clearvars;

% This script will test how the speed of the models grows with N0. the data
% generated here can be used to reproduce Figure 3FG as well as some SI
% figures

base_name = 'Figure3FG_';
n = 2; % number of initial tumor sizes to test
nsamps = 1; % number of samples to run for each method at each initial tumor size
save_output = false; % set to true if you want to save the output

%% parameters that are constant for all patients
N0 = round(logspace(log10(1),log10(10),n)); % initial number of tumor cells

%% dosing regimens
aFGFR3_start_min = 1; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1; % days between successive doses of anti-FGFR3
DoW_start = 6; % day of week to start on; 0 = monday (so 6 means the first simmed day is considered a Sunday and the first dose should be given on Day 1 (see aFGFR3_start_simday below)
censor_date = 2; % day in simulation to end

n_doses_aFGFR3 = 5; % max number of doses of FGFR3 to give


%% set file names
filename = nextFileName('data',base_name,3);

%% initialize pars
pars = basePars_FGFR3();
pars.main.max_pde_dt = (.25/60)/24;


aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min; % simulation day to start therapy
events = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]); % decide on when events occur within each simulation; events include new doses and censoring

f = @(ni,pars) startPatient_FGFR3(N0(ni),pars,events{1}); % function to run simulation with a given set of parameters

method = {'Local','Global'};
sz = [length(method),nsamps,n];

% preparing the timing of the runs
total_runs = prod(sz);
times = zeros(sz);

mu_n = 0; % weighted average of simulation times for displaying estimated time remaining for cohort
timer_start = tic;

for i = total_runs:-1:1
    [mi,si,ni] = ind2sub(sz,i); % this way it starts with the largest N0 samples and so will initially vastly overestimate the time remaining
    q = pars;
    q.main.method = method{mi};

    timerVal = tic;
    tracked(mi,si,ni) = f(ni,q); % run (and time) the simulation
    times(mi,si,ni) = toc(timerVal);

    % displays estimated time remaining for whole cohort (not part of the
    % method, but helpful to let you know how long to expect)
    num_done = total_runs - i + 1;
    mu_n = ((num_done-1)*mu_n+2*times(mi,si,ni))/(num_done+1); % this is computing the average duration of each run, weighting the more recent runs more heavily
    etr = mu_n * (total_runs-num_done);
    fprintf('Finished %d of %d, or %3.2f%%, after %s. ETR: %s for total run time of %s.\n',...
        num_done,total_runs,100*num_done/total_runs,duration(0,0,times(mi,si,ni)),duration(0,0,etr),duration(0,0,etr+toc(timer_start)))

end

if save_output
    if ~exist('./data','dir') % make sure the data directory exists
        mkdir('./data')
    end
    filename = nextFileName('data',base_name,3);
    save(filename,'-v7.3')
end

%% generate panels

nfigs = 2;
figs = gobjects(nfigs,1);
method_colors = lines(2);
lsty = {'-.','-'}; % {local,global} method line styles
t = tracked(1).T;
nt = length(t);

%% total wall time
figs(1) = figure("Name","Wall Time Dependence on Tumor Size");
hold on;
for mi = 1:2
    plot(N0,squeeze(mean(times(mi,:,:),1)),'Color',method_colors(mi,:),'LineWidth',2,'LineStyle',lsty{mi},'Marker','o','MarkerSize',20,'MarkerFaceColor',method_colors(mi,:))
end
set(gca,'YScale','log','FontSize',20)
legend({'Local','Global'},"Location","best")
xlabel('Initial Tumor Size')
ylabel('Wall Time(s)')

%% ode wall time proportion
total_wall_times = zeros([nt-1,sz]); % wall times for each update (no update made at last time point since the simulation ends there)
ode_wall_times = zeros([nt-1,sz]); % wall times for solving ODEs in each update (no update made at last time point since the simulation ends there)
for ni = 1:n
    for si = 1:nsamps
        for mi = 1:2
            total_wall_times(:,mi,si,ni) = tracked(mi,si,ni).wall_time;
            ode_wall_times(:,mi,si,ni) = tracked(mi,si,ni).wall_time_ode;
        end
    end
end
ode_start_ind = find(ode_wall_times(:,1,1,1)>0,1); % the first update step that required solving an ODE (after drug administered)

figs(2) = figure("Name","ODE Wall Time Proportion");
hold on;
for mi = 1:2
    plot(N0,squeeze(mean(ode_wall_times(ode_start_ind:end,mi,:,:),[1,3])),'Color',method_colors(mi,:),'LineWidth',2,'LineStyle',lsty{mi},'Marker','o','MarkerSize',20,'MarkerFaceColor',method_colors(mi,:));
end
set(gca,'YScale','log','FontSize',20)
legend({'Local','Global'},"Location","best")
xlabel('Initial Tumor Size')
ylabel('Wall Time(s)')


