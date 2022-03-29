function [tumors,aFGFR3] = updateTumor_FGFR3(p,order,tumors,V_tot,NT,...
    aFGFR3,location_inds,event_ind,proliferation_timer_ind,tumor_receptors_inds)

for ord_ind=1:length(order)
    j = order(ord_ind);
    switch tumors(j,event_ind)
        case 1 % proliferation
            locs=tumors(j,location_inds) + p.neighbors;
            trouble_tum = max(abs(tumors(:,location_inds)-tumors(j,location_inds)),[],2)==1;
            free_locs = locs(~any(all(locs == reshape(tumors(trouble_tum,location_inds)',1,3,[]),2),3),:); % see test_setdiff.m in Step22 for why this is the way I'm determining this
            if size(free_locs,1)>=(26-p.occmax) % check that the cell met the qualifications for proliferating
                NT = NT+1;
                distances = sum(abs(tumors(j,location_inds)-free_locs),2);
                ind = randsample(size(free_locs,1),1,true,1./distances);
                tumors(end+1,location_inds) = free_locs(ind,:); % put new tumor cell at random open location
                tumors([j,end],proliferation_timer_ind) = p.min_prolif_wait - (p.dt-tumors(j,proliferation_timer_ind))*rand(); % must wait min_prolif_wait days before proliferating; assume the proliferation happened sometime in the interval [tumors(j,proliferation_timer_ind),dt] (uniformly) so some progress towards next prolif has happened
                tumors([j,end],tumor_receptors_inds(1)) = tumors(j,tumor_receptors_inds)*[1;1;0;.5;1;1]; % monomer concentrations stay the same; all complexes get split evenly between the two, so to compensate half the number of monomers in them is added back to the monomer species
                tumors([j,end],tumor_receptors_inds([2,4:end])) = .5*tumors([j,j],tumor_receptors_inds([2,4:end])); % complexes end up split evenly between the two
                tumors(end,tumor_receptors_inds(3)) = aFGFR3; % set inhibitor concentration at new tumor to be ambient inhibitor concentration (this only matters for the global method; in local, aFGFR3 here is 0 so that it can be set just before next ODE step)
            else
                tumors(j,proliferation_timer_ind) = 0; % if not enough empty space, then allow this cell to try proliferating again
                tumors(j,event_ind) = 6; % say this cell experienced contact inhibition
            end % end of if
            
        case 2 % spontaneous apoptosis
            aFGFR3 = ((V_tot-NT)*aFGFR3+tumors(j,tumor_receptors_inds(3)))/(V_tot-NT+1); % take weighted average of ambient inhibitor concentration with the cleared one. Only matters for global method
            NT = NT-1;
            
        otherwise
            error('should not do nothing')
            
    end % end of switch
end % end of j for
