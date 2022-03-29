function [y,ambient] = fullGlobalODE_FGFR3(ic,p,dt,aFGFR3_circ0)

options = odeset('Nonnegative',1:numel(ic));

[~,Y] = ode15s(@(t,x) globalODE(t,x,p,aFGFR3_circ0),linspace(0,dt,3),ic(:),options);
y = Y(end,3:end-1);
ambient = Y(end,end);

function dx = globalODE(t,x,p,aFGFR3_circ0)

% x = [T;Y; [ RF,;DA;C;R_F_C;D_C_C;D_A_C;ambient drug ] ] molecules all given by average concentrations

dx = zeros(size(x));

% free FGFR3
dx(3) = -2*p.kf*x(3)^2 + 2*p.kr*x(4) + 2*p.kp*x(4) ... % dimerization and backwards reaction
    - p.k_on_R*x(3)*x(5) + p.k_off_R*x(6); % FGFR3 reacting with aFGFR3


% dimerized FGFR3
dx(4) = p.kf*x(3)^2 - p.kr*x(4) - p.kp*x(4) ... % dimerization and backwards reaction
    - p.k_on_D*x(4)*x(5) + p.k_off_D*x(8); % aFGFR3 reacting with dimers

% free aFGFR3
dx(5) = p.prop_tov_receiving_inhibitor*(p.aFGFR3_influx*aFGFR3_circ0*exp(-p.aFGFR3_sysdecay*t)-p.aFGFR3_eflux*x(5)) ... % aFGFR3 in circulation flowing into tumor-occupied volume and flowing back into circulation
    - p.aFGFR3_degradation*x(5) ... % aFGFR3 leaving the tumor-occupied volume
    - p.k_on_R*x(3)*x(5) + p.k_off_R*x(6) ...% aFGFR3-FGFR3-dimer reactions
    - p.k_on_D*x(4)*x(5) + p.k_off_D*x(8) ... % aFGFR3-FGFR3-dimer reactions
    + 2*(1-p.prop_tum_neighbors_that_are_tumors)*p.aFGFR3_diffusion*(x(9)-x(5))*(3/p.cell_width^2); % drug moving within TME

% aFGFR3-FGFR3 complex
dx(6) = p.k_on_R*x(3)*x(5) - p.k_off_R*x(6)... % FGFR3 reacting with aFGFR3
    - 2*p.kf*x(6)^2 + 2*p.kr*x(7) + 2*p.kp*x(7); % dimerization of these complexes (and the backwards reaction)

% aFGFR3-FGFR3 dimer
dx(7) = p.kf*x(6)^2 - p.kr*x(7) - p.kp*x(7); % dimerization of aFGFR3-FGFR3 complexes (and the backwards reaction)

% aFGFR3-dimer complex
dx(8) = p.k_on_D*x(4)*x(5) - p.k_off_D*x(8); % aFGFR3 reacting with dimers

% ambient aFGFR3 concentration
dx(9) = p.prop_amb_receiving_inhibitor*(p.aFGFR3_influx*aFGFR3_circ0*exp(-p.aFGFR3_sysdecay*t)-p.aFGFR3_eflux*x(9)) ... % influx from circulation and flowing back into circulation
    - p.aFGFR3_degradation*x(9) ... % degradation and eflux
    + 2*(1-p.prop_free_neighbors_that_are_free)*p.aFGFR3_diffusion*(x(5)-x(9))*(3/p.cell_width^2); % diffusion of drug within TME
