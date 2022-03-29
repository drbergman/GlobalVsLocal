function tracked = saveSubstrateInfo(tracked,solver,method,substrate,tumors_ind_ind,counter)

for si = 1:solver.num_substrates
    if solver.is_present(si)
        if method=="local"
            ambient_temp = substrate(si).me_concentration;
            ambient_temp(tumors_ind_ind) = [];

            tracked.substrate_ambient(counter,si) = mean(ambient_temp);
            tracked.substrate_TOV(counter,si) = mean(substrate(si).me_concentration(tumors_ind_ind));

            tracked.substrate_by_z(counter,:,si) = squeeze(mean(substrate(si).me_concentration,1:2));
            tracked.substrate_by_z_std(counter,:,si) = squeeze(std(substrate(si).me_concentration,[],1:2));
        elseif method=="global"
            conc_temp = substrate(si).me_concentration(solver.regions);

            tracked.substrate_TOV(counter,si) = mean(conc_temp(tumors_ind_ind));
            conc_temp(tumors_ind_ind) = [];

            tracked.substrate_ambient(counter,si) = mean(conc_temp);

            if solver.substrate_entry=="floor"
                tracked.substrate_by_z(counter,:,si) = substrate(si).me_concentration;
                tracked.substrate_by_z_std(counter,:,si) = 0;
            else
                tracked.substrate_by_z(counter,:,si) = reshape(mean(substrate(si).me_concentration(solver.regions),1:2),1,[],1);
                tracked.substrate_by_z_std(counter,:,si) = reshape(std(substrate(si).me_concentration(solver.regions),[],1:2),1,[],1);
            end
        else
            error("unknown method")
        end
    else
        tracked.substrate_ambient(counter,si) = 0;
        tracked.substrate_TOV(counter,si) = 0;
        if solver.is_pk(si)
            tracked.substrate_circ(counter,si) = 0;
            tracked.substrate_periph(counter,si) = 0;
        end
        tracked.substrate_by_z(counter,:,si) = 0;
        tracked.substrate_by_z_std(counter,:,si) = 0;
    end
end