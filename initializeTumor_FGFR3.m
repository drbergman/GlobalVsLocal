function tumors = initializeTumor_FGFR3(N0,RT,location_inds,...
    event_ind,tumor_clearance_inds,tumor_receptors_inds,proliferation_timer_ind)

radius = nthroot((3/4)*(1/pi())*N0,3); % radius of sphere with volume equal to the number of cells
s = ceil(4*radius); % size of TME to hold the simulation
side = floor((1-s/2):(s/2));
sidebar = mean(side);
locs = allCombos(side,side,side);

tumors = zeros(N0,tumor_receptors_inds(end));

tumors(:,location_inds) = locs(datasample(1:s^3,N0,'Replace',false,'Weights',exp(-1*sqrt(  sum((locs-sidebar).^2,2)  ))),:);
tumors(:,tumor_clearance_inds(1)) = 0; % proportion cleared by immune cell (if it reaches 1, then the tumor is marked for apoptosis; if its clearance rate is set to 0, then it is cleared with this probability)
tumors(:,tumor_clearance_inds(2)) = 0; % rate at which the tumor is being cleared by the immune system.
tumors(:,event_ind) = 0; % event index
tumors(:,proliferation_timer_ind) = 0; % time until next possible proliferation (days)

tumors(:,tumor_receptors_inds(1)) = RT; % inactive monomers of FGFR3 (nM)
tumors(:,tumor_receptors_inds(2)) = 0; % active dimers of signaling FGFR3 (nM)
tumors(:,tumor_receptors_inds(3)) = 0; % inhibitor (nM)
tumors(:,tumor_receptors_inds(4)) = 0; % inhibitor-monomer complex (nM)
tumors(:,tumor_receptors_inds(5)) = 0; % inhibitor-monomer dimer (nM)
tumors(:,tumor_receptors_inds(6)) = 0; % inhibitor-dimer complex (nM)