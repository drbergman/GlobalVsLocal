function dx = agentODE_FGFR3(x,kf,kr,kp,k_on_R,k_off_R,k_on_D,k_off_D)

dx = zeros(size(x));

% free FGFR3
dx(:,1) = -2*kf*x(:,1).^2 + 2*(kr+kp)*x(:,2)...% dimerization and backwards reaction
    - k_on_R*x(:,1).*x(:,3) + k_off_R*x(:,4); % FGFR3 reacting with aFGFR3

% dimerized FGFR3
dx(:,2) = kf*x(:,1).^2 - (kr+kp)*x(:,2) ...% dimerization and backwards reaction
    - k_on_D*x(:,2).*x(:,3) + k_off_D*x(:,6);% aFGFR3 reacting with dimers

% free aFGFR3
dx(:,3) = - k_on_R*x(:,1).*x(:,3) + k_off_R*x(:,4)...% aFGFR3-FGFR3-dimer reactions
    - k_on_D*x(:,2).*x(:,3) + k_off_D*x(:,6); % aFGFR3-FGFR3-dimer reactions

% aFGFR3-FGFR3 complex
dx(:,4) = k_on_R*x(:,1).*x(:,3) - k_off_R*x(:,4)...% FGFR3 reacting with aFGFR3
    - 2*kf*x(:,4).^2 + 2*(kr+kp)*x(:,5); % dimerization of these complexes (and the backwards reaction)

% aFGFR3-FGFR3 dimer
dx(:,5) = kf*x(:,4).^2 - (kr+kp)*x(:,5); % dimerization of aFGFR3-FGFR3 complexes (and the backwards reaction)

% aFGFR3-dimer complex
dx(:,6) = k_on_D*x(:,2).*x(:,3) - k_off_D*x(:,6); % aFGFR3 reacting with dimers