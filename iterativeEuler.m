function EC = iterativeEuler(IC,t,dt,fn,p)

% Euler method with nonnegativity constraint: if update results in negative
% state variable, subdivide interval in half and repeat

dx = fn(t,IC,p);
EC = IC + dt*dx;

if any(EC<0)
    EC1 = iterativeEuler(IC,t,0.5*dt,fn,p);
    EC = iterativeEuler(EC1,t+0.5*dt,0.5*dt,fn,p);
end
