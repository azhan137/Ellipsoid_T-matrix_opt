%take the first derivative of bessel function of first kind or hankel of
%the first kind
%nu = 1: bessel
%nu = 3: hankel
%Z: bessel argument
%l: bessel order


%d(j_(l))/dz = (l+1)j_(l)-z*j_(l+1)

function first_derivative = d1Z_sph_bessel(nu,l,Z)

first_derivative = l./Z.*sph_bessel(nu,l,Z)-sph_bessel(nu,l+1,Z);