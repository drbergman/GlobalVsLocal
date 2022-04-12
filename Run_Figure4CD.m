clearvars;

% This script will create the drug concentration spheres from Figure 4CD.


base_name = 'Figure4CD_';
nsamps = 1; % number of samples to run for each method
N0 = 1e4; % initial number of tumor cells
save_output = false; % set to true if you want to save the output

%% dosing regimens
aFGFR3_start_min = 0; % min day in simulation to start. will pick next Monday on or after this to start
aFGFR3_days_between = 1; % days between successive doses of anti-FGFR3
DoW_start = 6; % day of week to start on; 0 = monday (so 6 means the first simmed day is considered a Sunday and the first dose should be given on Day 1 (see aFGFR3_start_simday below)
censor_date = 8; % day in simulation to end

n_doses_aFGFR3 = 5; % max number of doses of FGFR3 to give

%% initialize pars
pars = basePars_FGFR3();
pars.plot_properties.plot_inhibitor = true; % make sure to plot the inhibitor
pars.main.method = 'Local'; % only running the local method for this

aFGFR3_start_simday = 7-(mod(aFGFR3_start_min + DoW_start-1,7)+1) + aFGFR3_start_min; % simulation day to start therapy

events = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,[],[]); % decide on when events occur within each simulation; events include new doses and censoring

startPatient_FGFR3(N0,pars,events{1}); % function to run simulation with a given set of parameters
