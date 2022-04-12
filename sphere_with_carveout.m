function [sx,sy,sz,sc] = sphere_with_carveout(x,y,z,c)

% sphere_with_carveout will plot a sphere with the first octant carved out
% so you can see into the sphere. 
% x: vector of x coordinates
% y: vector of y coordinates
% z: vector of z coordinates
% c: color data on ndgrid(x,y,z)

[sx,sy,sz] = sphere(200);
first_octant_indices = sx>=0 & sy>=0 & sz>=0;

% [X,Y,Z] = meshgrid(linspace(-1,1),linspace(-1,1),linspace(-1,1));
% omega = [2,3,4];
% TME = cos(omega(1)*pi*X).*cos(omega(2)*pi*Y).*cos(omega(3)*pi*Z);
% TME = sqrt(X.^2+Y.^2+Z.^2);

[sx(first_octant_indices),sy(first_octant_indices),sz(first_octant_indices)] = ...
    carveout(sx(first_octant_indices),sy(first_octant_indices),sz(first_octant_indices));

sx = (x(end)-x(1))*.5*sx + mean(x); % make the sphere the appropriate size
sy = (y(end)-y(1))*.5*sy + mean(y); % make the sphere the appropriate size
sz = (z(end)-z(1))*.5*sz + mean(z); % make the sphere the appropriate size
sc = interp3(y,x,z,c,sx,sy,sz,'spline');
end
function [x,y,z] = carveout(x,y,z)

R{1} = z>y & x>=y;
R{2} = y>x & z>x;
R{3} = x>=z & y>=z;

for i = 1:3
    xR{i} = x(R{i});
    yR{i} = y(R{i});
    zR{i} = z(R{i});
end

m = [1;1;1]/sqrt(3);
v0 = [[1;0;0],[0;1;0],[1;0;0]]; % a vector on each plane
% v = {[1;0;0]-[[0;0;1],m],...
%     [0;0;1]-[[0;1;0],m],...
%     [1;0;0]-[[0;1;0],m]}; % vectors in direction of each plane I want to project onto

% v = {[[0;0;1],[1;0;0]]-m,...
%     [[0;1;0],[0;0;1]]-m,...
%     [[1;0;0],[0;1;0]]-m}; % vectors in direction of each plane I want to project onto
M{1} = [[0;0;1],[1;0;0],m];
M{2} = [[0;1;0],[0;0;1],m];
M{3} = [[1;0;0],[0;1;0],m];

u = [1/sqrt(2) * sqrt((1-1/sqrt(3))/2);1/sqrt(2) * sqrt((1-1/sqrt(3))/2);sqrt((1+1/sqrt(3))/2)]; % axis of rotation to bring m to [0 0 1]
Rm = -eye(3) + 2*(u*u'); % rotation matrix to move m to [0 0 1]

ex_phi = linspace(0,pi/2,1000);
ex_r = acos(Rm(3,:)*[cos(ex_phi);sin(ex_phi);zeros(1,length(ex_phi))]);
for i = 1:3
%     P{i} = v{i}*((v{i}'*v{i})\(v{i}'*([xR{i},yR{i},zR{i}]'-v0(:,i))))+v0(:,i); % the projection onto the three planes
%     BC{i} = M{i}\P{i}; % Barycentric coordinates on these planes
%     NC{i} = M{i}(:,1:2)*BC{i}(1:2,:); % new coordinates
%     x(R{i}) = NC{i}(1,:);
%     y(R{i}) = NC{i}(2,:);
%     z(R{i}) = NC{i}(3,:);
    
    rotated_points{i} = Rm*[xR{i},yR{i},zR{i}]';
    r{i} = acos(rotated_points{i}(3,:));
    p0 = Rm*M{i}(:,1);
    pf = Rm*M{i}(:,2);
    phi0{i} = cart2pol(p0(1)/sqrt(1-p0(3)^2),p0(2)/sqrt(1-p0(3)^2));
    phif{i} = cart2pol(pf(1)/sqrt(1-pf(3)^2),pf(2)/sqrt(1-pf(3)^2));
    if phif{i}<phi0{i}
        phif{i} = phif{i}+2*pi;
    end
%     phi0{i} = acos(p0(1)/sqrt(1-p0(3)^2));
%     if abs(sin(phi0{i})-p0(2)/sqrt(1-p0(3)^2))>1e-5
%         phi0{i} = 2*pi-phi0{i};
%     end
%     phif{i} = acos(pf(1)/sqrt(1-pf(3)^2));
%     if abs(sin(phif{i})-pf(2)/sqrt(1-pf(3)^2))>1e-5
%         phif{i} = 2*pi-phif{i};
%     end
%     phi{i} = (pi/2)*(acos(rotated_points{i}(1,:)./sin(r{i})) - phi0{i})./(phif{i}-phi0{i});
    base_angle = cart2pol(rotated_points{i}(1,:)./sin(r{i}),rotated_points{i}(2,:)./sin(r{i}));
    base_angle(base_angle<phi0{i}) = base_angle(base_angle<phi0{i}) + 2*pi;
    phi{i} = (pi/2)*(base_angle-phi0{i})./(phif{i}-phi0{i});
%     NC{i} = interp1([min(r{i});max(r{i})],[0;1],r{i}).*(M{i}(:,1:2)*[cos(phi{i});sin(phi{i})]);
    NC{i} = r{i}./interp1(ex_phi,ex_r,phi{i}).*(M{i}(:,1:2)*[cos(phi{i});sin(phi{i})]);
    x(R{i}) = NC{i}(1,:);
    y(R{i}) = NC{i}(2,:);
    z(R{i}) = NC{i}(3,:);
end
end
    
% end
% 
% figure; hold on
% for i = 1:3
% %     scatter3(P{i}(1,:),P{i}(2,:),P{i}(3,:),'filled')
%     scatter3(NC{i}(1,:),NC{i}(2,:),NC{i}(3,:),'filled')
% end
% xlabel('x')
% ylabel('y')
% zlabel('z')