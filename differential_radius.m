





%for an elliosoid, returns gradient of ellipsoid surface with respect to
%axes a (1), b (2), c (3) at a given angular position (theta, phi)
%rotation about x,y plane (a,b) is defined as phi (4)

function grad_r = differential_radius(theta,phi,a,b,c)

st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);

r = 1./sqrt(st.^2.*(cp.^2./a^2+sp.^2./b^2)+ct.^2./c^2);

grad_r = zeros([size(theta),4]);

prefactor = r.^3.*st.^2;

grad_r(:,:,1) = prefactor.*cp.^2/a^3;
grad_r(:,:,2) = prefactor.*sp.^2/b^3;
grad_r(:,:,3) = r.^3.*ct.^2/c^3;
grad_r(:,:,4) = -1*prefactor.*sp.*cp.*(1/b^2-1/a^2);
