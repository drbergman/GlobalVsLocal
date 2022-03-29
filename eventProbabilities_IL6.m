function p_matrix = eventProbabilities_IL6(tumors,event_pars,dt,inds,phiD)

r_matrix = [event_pars.prolif_rate(tumors(:,inds.type_ind))',... % proliferation rates
    (event_pars.delta(tumors(:,inds.type_ind))')./(1+phiD/event_pars.ec50),... % spontaneous apoptosis rates
    event_pars.move_rate + zeros(size(tumors,1),1)]; % movement rates

p_matrix = 1-exp(-r_matrix.*[max(0,dt-tumors(:,inds.proliferation_timer_ind)),...
                            dt*ones(size(tumors,1),size(r_matrix,2)-1)]); % use dt as time step except for proliferation, use dt - remaining time to wait before next proliferation

assert(all(sum(p_matrix,2)<=1)) % make sure the probabilities do not add up to more than 1
