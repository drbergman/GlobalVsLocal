function Js = localJacobian_vectorized(Js,x,kf,kr,kp,k_on_R,k_off_R,k_on_D)

Js(1,1,:) = -4*kf*x(:,1)-k_on_R*x(:,3);
Js(1,3,:) = -k_on_R*x(:,1);

Js(2,1,:) = 2*kf*x(:,1);
Js(2,2,:) = -kr-kp-k_on_D*x(:,3);
Js(2,3,:) = -k_on_D*x(:,2);

Js(3,1,:) = -k_on_R*x(:,3);
Js(3,2,:) = -k_on_D*x(:,3);
Js(3,3,:) = -k_on_R*x(:,1)-k_on_D*x(:,2);

Js(4,1,:) = k_on_R*x(:,3);
Js(4,3,:) = k_on_R*x(:,1);
Js(4,4,:) = -k_off_R-4*kf*x(:,4);

Js(5,4,:) = 2*kf*x(:,4);

Js(6,2,:) = k_on_D*x(:,3);
Js(6,3,:) = k_on_D*x(:,2);