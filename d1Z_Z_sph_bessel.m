%take the first derivative of bessel function of first kind or hankel of
%the first kind
%nu = 1: bessel
%nu = 3: hankel
%Z: bessel argument
%l: bessel order


%d(z*j_(l))/dz = (l+1)j_(l)-z*j_(l+1)

function first_derivative = d1Z_Z_sph_bessel(nu,l,Z)
if l ~= 1
    first_derivative = (l+1)*sph_bessel(nu,l,Z)-Z.*sph_bessel(nu,l+1,Z);
else
    first_derivative = Z.*sph_bessel(nu,0,Z)-sph_bessel(nu,1,Z);
end