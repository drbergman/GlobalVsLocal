function [tumors,L] = updateTumor_IL6(tumors,L,substrate,order,inds,grid,update_pars,method,phiD,solver)

NT = size(tumors,1);
for ord_ind = 1:length(order)
    j = order(ord_ind); % index of current cell in tumors

    switch(tumors(j,inds.event_ind))

        case 1 % proliferation
            n_ind = tumors(j,inds.ind_ind)+grid.rel_pos_ind; % neighbor indices
            if (any(tumors(j,inds.subs_inds)==1) || any(tumors(j,inds.subs_inds)==grid.size))

                % enforce reflecting boundary conditions

                % if on border in x direction
                if tumors(j,inds.subs_inds(1))==1
                    n_ind([1:3:13,15:3:24]) = tumors(j,inds.ind_ind);
                end
                if tumors(j,inds.subs_inds(1))==grid.size(1)
                    n_ind([3:3:12,14:3:26]) = tumors(j,inds.ind_ind);
                end

                % if on border in y direction
                if tumors(j,inds.subs_inds(2))==1
                    n_ind([1:3,10:12,18:20]) = tumors(j,inds.ind_ind);
                end
                if tumors(j,inds.subs_inds(2))==grid.size(2)
                    n_ind([7:9,15:17,24:26]) = tumors(j,inds.ind_ind);
                end

                % if on border in z direction
                if tumors(j,inds.subs_inds(3))==1
                    n_ind(1:9) = tumors(j,inds.ind_ind);
                end
                if tumors(j,inds.subs_inds(3))==grid.size(3)
                    n_ind(18:26) = tumors(j,inds.ind_ind);
                end

            end

            if nnz(L(n_ind))<=update_pars.occmax % check how many neighbors are occupied
                ind = randsample(26,1,true,(~L(n_ind))./sqrt(sum(update_pars.neighbors.^2,2))); % weight by whether sites are empty and by the reciprocal of the distance
                rel_loc = update_pars.neighbors(ind,:); % get the relative position of the new spot
                tumors(end+1,inds.location_inds) = tumors(j,inds.location_inds)+rel_loc;
                
                tumors(end,inds.subs_inds) = tumors(j,inds.subs_inds)+rel_loc;
                assert(tumors(end,inds.subs_inds(3))>=1)
                
                tumors(end,inds.ind_ind) = n_ind(ind);
                L(n_ind(ind)) = L(tumors(j,inds.ind_ind)); % set value at lattice site based on original cell (this should work instead of above commented out line)
                tumors([j,end],inds.proliferation_timer_ind) = update_pars.min_prolif_wait - (update_pars.dt-tumors(j,inds.proliferation_timer_ind))*rand();
                if method=="local"
                    for si = 1:2
                        tumors(end,inds.tumor_receptors_inds(solver.substrate_inds(si))) = substrate(si).me_concentration(n_ind(ind));
                    end
                elseif method=="global"
                    for si = 1:2
                        tumors(end,inds.tumor_receptors_inds(solver.substrate_inds(si))) = substrate(si).me_concentration(solver.regions(n_ind(ind)));
                    end
                else
                    error("unknown method")
                end

                if tumors(j,inds.type_ind)==2
                    if tumors(j,inds.num_prolifs) == 1 % then the progenitor cell divides into two terminally differentiated cells
                        tumors([j,end],inds.type_ind) = 3;
                    else
                        tumors(end,inds.type_ind) = 2;
                        tumors([j,end],inds.num_prolifs) = tumors(j,inds.num_prolifs) - 1;
                    end
                elseif tumors(j,inds.type_ind)==1
                    psmin = update_pars.mu_S * (update_pars.P_Smax - update_pars.P_Smin_star)*phiD(j) + update_pars.P_Smin_star;
                    psymm = (update_pars.P_Smax - psmin) * update_pars.P_ec50^update_pars.hill_coeff / (update_pars.P_ec50^update_pars.hill_coeff+sum(tumors(:,inds.type_ind)==1)^update_pars.hill_coeff) + psmin;
                    if rand() < psymm
                        tumors(end,inds.type_ind) = 1;
                    else % then asymmetric division
                        tumors(end,inds.type_ind) = 2;
                        tumors(end,inds.num_prolifs) = update_pars.w;
                        tumors(end,inds.tumor_receptors_inds(1)) = update_pars.RT(2);
                        continue; % don't need to set bound complex concentrations
                    % P: probability of stem cell self-renewal
                    % P_Smin: minimum probability of stem cell self-renewal
                    % S: number of stem cells
                    % phi_S: fractional occupancy of stem cells
                    % P_Smin = @(phi_S) mu_S * (P_Smax - P_Smin_star)*phi_S + P_Smin_star;
                    % P = @(S,phi_S) = (P_Smax - P_Smin(phi_S)) * P_ec50^hill_coeff / (P_ec50^hill_coeff+S^hill_coeff) + P_Smin(phi_S);
                    end
                end
                tumors(j,inds.tumor_receptors_inds([3,5])) = 0.5*tumors(j,inds.tumor_receptors_inds([3,5])); % split complexes evenly between two daughter cells
                tumors(end,inds.tumor_receptors_inds([3,5])) = tumors(j,inds.tumor_receptors_inds([3,5])); % split complexes evenly between two daughter cells
                tumors([j,end],inds.tumor_receptors_inds(1)) = tumors(j,inds.tumor_receptors_inds(1)) + sum(tumors(j,inds.tumor_receptors_inds([3,5]))); % the IL6R lost above become unbound and split evenly between the two
                
            else
                tumors(j,inds.proliferation_timer_ind) = 0; % if not enough empty space, then allow this cell to try proliferating again
            end

        case 2 % apoptosis
            NT = NT-1;

        case 3 % movement
            n_ind = tumors(j,inds.ind_ind)+grid.rel_pos_ind; % neighbor indices

            if (any(tumors(j,inds.subs_inds)==1) || any(tumors(j,inds.subs_inds)==grid.size))
                % directions of proliferating are limited; set any
                % direction that moves off grid to center point. this way,
                % those spots are disallowed from proliferating into

                % enforce reflecting boundary conditions
                
                % if on border in x direction
                if tumors(j,inds.subs_inds(1))==1
                    n_ind([1:3:13,15:3:24]) = tumors(j,inds.ind_ind);
                end
                if tumors(j,inds.subs_inds(1))==grid.size(1)
                    n_ind([3:3:12,14:3:26]) = tumors(j,inds.ind_ind);
                end

                % if on border in y direction
                if tumors(j,inds.subs_inds(2))==1
                    n_ind([1:3,10:12,18:20]) = tumors(j,inds.ind_ind);
                end
                if tumors(j,inds.subs_inds(2))==grid.size(2)
                    n_ind([7:9,15:17,24:26]) = tumors(j,inds.ind_ind);
                end

                % if on border in z direction
                if tumors(j,inds.subs_inds(3))==1
                    n_ind(1:9) = tumors(j,inds.ind_ind);
                end
                if tumors(j,inds.subs_inds(3))==grid.size(3)
                    n_ind(18:26) = tumors(j,inds.ind_ind);
                end

            end

            if nnz(L(n_ind))<26 % make sure there is a direction to move
                L(tumors(j,inds.ind_ind)) = false; % remove cell from lattice (if a direction would have moved the cell off the lattice, then that selection will now move it here, i.e. no movement)
                
                ind = randsample(26,1,true,(~L(n_ind))./sqrt(sum(update_pars.neighbors.^2,2))); % weight by whether sites are empty and by the reciprocal of the distance
                
                if n_ind(ind)==tumors(j,inds.ind_ind)
                    % cell ended up staying in place
                    L(tumors(j,inds.ind_ind)) = true;
                    continue;
                end
                rel_loc = update_pars.neighbors(ind,:); % get the relative position of the new spot
                tumors(j,inds.location_inds) = tumors(j,inds.location_inds)+rel_loc;
                
                assert(tumors(j,inds.subs_inds(3))>=1)
                tumors(j,inds.subs_inds) = tumors(j,inds.subs_inds)+rel_loc;
                assert(tumors(j,inds.subs_inds(3))>=1)
                
                tumors(j,inds.ind_ind) = n_ind(ind);
                L(n_ind(ind)) = true; % put cell back on lattice
                if method=="local"
                    tumors(j,inds.tumor_receptors_inds(2)) = substrate(1).me_concentration(n_ind(ind)); % split complexes evenly between two daughter cells
                    tumors(j,inds.tumor_receptors_inds(4)) = substrate(2).me_concentration(n_ind(ind)); % split complexes evenly between two daughter cells
                elseif method=="global"
                    for si = 1:2
                        tumors(j,inds.tumor_receptors_inds(solver.substrate_inds(si))) = substrate(si).me_concentration(solver.regions(n_ind(ind)));
                    end
                else
                    error("unknown method")
                end
            end

        otherwise
            error('picked an event that I have not yet coded')
    end

end