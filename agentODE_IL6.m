function dx = agentODE_IL6(~,x,p)

dx = (-p.kf*x(1,:)*x(2,:)+p.kr*x(3,:)).*[1;1;-1;0;0] + (-p.kf_I*x(1,:)*x(4,:) + p.kr_I*x(5,:)).*[1;0;0;1;-1];
dx(2,:) = dx(2,:) + p.secretion;