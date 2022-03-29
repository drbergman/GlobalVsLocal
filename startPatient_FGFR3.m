function tracked = startPatient_FGFR3(tumors,q,E)

%% initialize inputs
tracked = [];
aFGFR3_tracker = {0,0};
grids = [];
current_day = 0;

DT = q.main.dt;

q.main.receptors_without_inhibitor = (q.main.kr+q.main.kp - sqrt((q.main.kr+q.main.kp)^2+8*q.main.kf*(q.main.kr+q.main.kp)*q.main.RT))/(-4*q.main.kf); % if no anti-fgfr3 therapy, then all tumor cells quickly reach this concentration of dimer
q.main.all_receptors_without_inhibitor = [q.main.receptors_without_inhibitor,(q.main.RT-q.main.receptors_without_inhibitor)/2,zeros(1,4)];   % if no anti-fgfr3 therapy, then all tumor cells quickly reach this concentration of fgfr3 stuff

for ei = 1:size(E,1)
    Nsteps = ceil((E(ei,1)-current_day)/DT);
    
    %% this block is for just in case dt needs to change to accomodate the event step
    q.main.dt = (E(ei,1)-current_day)/Nsteps;
    
    %% simulate
    [tracked,tumors,aFGFR3_tracker,grids] = ...
        simTumor_FGFR3(tumors,q,Nsteps,aFGFR3_tracker,grids,tracked);
    
    current_day = E(ei,1);
    
    %% update simulation based on next event
    switch E(ei,2)
        
        case 2
            aFGFR3_tracker{1} = aFGFR3_tracker{1}+q.main.aFGFR3_circ0;
            
        case Inf
            if ei<size(E,1)
                warning('censored before dosing all done')
            end
            break;
            
    end
    
end
end