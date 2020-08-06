%compute J matrices 11,12,21,22 using gaussian quadrature for computation
%   of full T matrix of an ellipsoid
%formulas are taken from 
%'Differential Cross Section of a Dielectric Ellipsoid by the T-Matrix 
%   Extended Boundary Condition Method', J.Schneider, and I. Peden. IEEE 
%   Transactions on Antennas and Propagation. Vol 36, NO 9, September 1988
% note that there are numerous typos in the above

%lmax: integer specifying maximum order
%Ntheta: number of points for Gaussian quadrature along polar
%Nphi: number of points for Gaussian quadrature along azimuthal
%a,b,c: axes of ellipsoid
%ni,ns: refractive index of medium and particle respectively
%lambda: wavelength of light
%nu: internal (1) or external (3)
%

function [J11,J12,J21,J22,dJ11,dJ12,dJ21,dJ22] = compute_J_ellip_celes_posm(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,nu)
%preallocate memory for J (nmax/2 x nmax/2)
nmax = jmult_max(1,lmax);
J11 = zeros(nmax/2);
J12 = zeros(nmax/2);
J21 = zeros(nmax/2);
J22 = zeros(nmax/2);
%preallocate for dJ (nmax/2 x nmax/2 x 3) 1:a, 2:b, 3:c
dJ11 = zeros(nmax/2,nmax/2,3);
dJ12 = zeros(nmax/2,nmax/2,3);
dJ21 = zeros(nmax/2,nmax/2,3);
dJ22 = zeros(nmax/2,nmax/2,3);
%setup for gaussian quadrature in two dimensions
[theta_i,wtheta_i] = generate_gauss_weights_abscissae(Ntheta,0,pi/2);
[phi_i,wphi_i] = generate_gauss_weights_abscissae(Nphi,0,pi);
%generate r, theta, phi corresponding to the surface of the ellisoid
[theta_map,phi_map] = meshgrid(theta_i,phi_i);
%r = ellip_rad(theta_map,phi_map,a,b,c);
%generate weight map
[wt,wp] = meshgrid(wtheta_i,wphi_i);
weight_map = wt.*wp;
%k vector freespace and in scatterer
k = 2*pi*ni/lambda;
ks = 2*pi*ns/lambda;
%associated legendre polynomials
P_lm = legendre_normalized_angular(theta_map,lmax);
[Pi_lm,Tau_lm] = spherical_functions_angular(theta_map,lmax);
ct = cos(theta_map);
st = sin(theta_map);
cp = cos(phi_map);
sp = sin(phi_map);
term = ((cp/a).^2+(sp/b).^2).*st.^2+(ct/c).^2;
r = 1./sqrt(term);
ellip1 = ellipsoid_function(a,b,c,phi_map,1).*st.*ct;
ellip2 = ellipsoid_function(a,b,c,phi_map,2).*sp.*cp.*st;
radius_grad = differential_radius(theta_map,phi_map,a,b,c);
[dellip1,dellip2] = ellipsoid_function_deriv(st,ct,sp,cp,a,b,c);

for li = 1:lmax
    %precompute bessel functions for li,mi and derivative
    b_li = sph_bessel(nu,li,k*r);
    db_li = d1Z_Z_sph_bessel(nu,li,k*r);
    db2_li = (li+li.^2-(k*r).^2)./(k*r).*b_li;
    d1b_li = d1Z_sph_bessel(nu,li,k*r);
    for mi = -li:li
        %compute index ni
        ni = multi2single_index(1,1,li,mi,lmax);
        %get spherical harmonics (p,pi,tau) and derivatives
        p_limi = P_lm{li+1,abs(mi)+1};
        pi_limi = Pi_lm{li+1,abs(mi)+1};
        tau_limi = Tau_lm{li+1,abs(mi)+1};
        
        for lp = 1:lmax
            %precompute bessel functions for lp,mp, and derivative
            j_lp = sph_bessel(1,lp,ks*r);
            dj_lp = d1Z_Z_sph_bessel(1,lp,ks*r);
            dj2_lp = (lp+lp.^2-(ks*r).^2)./(ks*r).*j_lp;
            d1j_lp = d1Z_sph_bessel(1,lp,ks*r);
            lfactor = 1/2/sqrt(li*(li+li)*lp*(lp+1));
            for mp = -lp:lp
                %compute index np
                np = multi2single_index(1,1,lp,mp,lmax);
                %get spherical harmonics (p,pi,tau) and derivatives
                p_lpmp = P_lm{lp+1,abs(mp)+1};
                pi_lpmp = Pi_lm{lp+1,abs(mp)+1};
                tau_lpmp = Tau_lm{lp+1,abs(mp)+1};
              
                %setup selection rules for J11,J22
                selection_rules_1122 = selection_rules(li,mi,lp,mp,1);
                
                %setup selection rules for J12,J21
                selection_rules_1221 = selection_rules(li,mi,lp,mp,2);
                
                %phi exponential phase factor
                phi_exp = exp(1i*(mp-mi)*phi_map);
                
                %compute J11,J22
                if selection_rules_1122 ~= 0
                    prefactor = selection_rules_1122*lfactor*phi_exp;
                    ang = mp*pi_lpmp.*tau_limi+mi*pi_limi.*tau_lpmp;
                    J11(ni,np) = sum(sum(-1i*weight_map.*r.^2.*st.*prefactor.*j_lp.*b_li.*(ang)));
                    dJ11_r = -1i*weight_map.*prefactor.*r.*st.*ang.*(r.*(k*d1b_li.*j_lp+ks*d1j_lp.*b_li)+2*b_li.*j_lp);
                    dJ11(ni,np,1) = sum(sum(dJ11_r.*radius_grad(:,:,1)));
                    dJ11(ni,np,2) = sum(sum(dJ11_r.*radius_grad(:,:,2)));
                    dJ11(ni,np,3) = sum(sum(dJ11_r.*radius_grad(:,:,3)));
                    J22_r = -1i*dj_lp.*db_li/ks/k.*st.*ang;
                    J22_t = 1i*r.^2/ks/k.*st.*(mp*li*(li+1)*dj_lp.*b_li.*p_limi.*pi_lpmp+mi*lp*(lp+1)*j_lp.*db_li.*p_lpmp.*pi_limi);
                    J22_p = r.^2/k/ks.*st.*(lp*(lp+1)*j_lp.*p_lpmp.*db_li.*tau_limi-li*(li+1)*dj_lp.*tau_lpmp.*b_li.*p_limi);
                    J22(ni,np) = sum(sum(weight_map.*prefactor.*(J22_r+J22_t.*ellip1+J22_p.*ellip2)));
                    dJ22_r = -1i/ks/k*st.*ang.*(k*db2_li.*dj_lp+ks*dj2_lp.*db_li);
                    dJ22_t = 1i*ellip1.*st/ks/k.*((lp+1)*lp*mi*p_lpmp.*pi_limi.*(2*r.*db_li.*j_lp+r.^2.*(k*db2_li.*j_lp+ks*d1j_lp.*db_li))+li*(li+1)*mp*p_limi.*pi_lpmp.*(2*r.*b_li.*dj_lp+r.^2.*(k*d1b_li.*dj_lp+ks*b_li.*dj2_lp)));
                    dJ22_p = ellip2.*st/ks/k.*(lp*(lp+1)*p_lpmp.*tau_limi.*(2*r.*db_li.*j_lp+r.^2.*(k*db2_li.*j_lp+ks*d1j_lp.*db_li))-li*(li+1)*p_limi.*tau_lpmp.*(2*r.*b_li.*dj_lp+r.^2.*(k*d1b_li.*dj_lp+ks*b_li.*dj2_lp)));
                    dJ22(ni,np,1) = sum(sum(weight_map.*prefactor.*((dJ22_r+dJ22_t+dJ22_p).*radius_grad(:,:,1)+J22_t.*dellip1(:,:,1)+J22_p.*dellip2(:,:,1))));
                    dJ22(ni,np,2) = sum(sum(weight_map.*prefactor.*((dJ22_r+dJ22_t+dJ22_p).*radius_grad(:,:,2)+J22_t.*dellip1(:,:,2)+J22_p.*dellip2(:,:,2))));
                    dJ22(ni,np,3) = sum(sum(weight_map.*prefactor.*((dJ22_r+dJ22_t+dJ22_p).*radius_grad(:,:,3)+J22_t.*dellip1(:,:,3))));
                end
                %compute J12,J21
                if selection_rules_1221 ~= 0
                    prefactor = selection_rules_1221*lfactor*phi_exp;
                    ang1 = mp*mi.*pi_lpmp.*pi_limi+tau_lpmp.*tau_limi;
                    term1 = ang1.*r.*st/k.*db_li.*j_lp;
                    term2 = -r.^3/k*li*(li+1).*st.*b_li.*p_limi.*j_lp.*(ellip1.*tau_lpmp+ellip2.*pi_lpmp*1i*mp);
                    J12(ni,np) = sum(sum(weight_map.*prefactor.*(term1+term2)));
                    dJ12_r = ang1/k.*st.*(j_lp.*db_li+r.*(ks*d1j_lp.*db_li+k*db2_li.*j_lp));
                    ang_t = -li*(li+1)*p_limi.*tau_lpmp/k.*st;
                    ang_p = -1i*li*(li+1)*mp*p_limi.*pi_lpmp/k.*st;
                    r3jb = r.^3.*j_lp.*b_li;
                    d_r3jb = 3*r.^2.*j_lp.*b_li+r.^3.*(k*d1b_li.*j_lp+ks*d1j_lp.*b_li);
                    dJ12(ni,np,1) = sum(sum(weight_map.*prefactor.*((dJ12_r+(ang_t.*ellip1+ang_p.*ellip2).*d_r3jb).*radius_grad(:,:,1)+r3jb.*(ang_t.*dellip1(:,:,1)+ang_p.*dellip2(:,:,1)))));
                    dJ12(ni,np,2) = sum(sum(weight_map.*prefactor.*((dJ12_r+(ang_t.*ellip1+ang_p.*ellip2).*d_r3jb).*radius_grad(:,:,2)+r3jb.*(ang_t.*dellip1(:,:,2)+ang_p.*dellip2(:,:,2)))));
                    dJ12(ni,np,3) = sum(sum(weight_map.*prefactor.*((dJ12_r+(ang_t.*ellip1+ang_p.*ellip2).*d_r3jb).*radius_grad(:,:,3)+r3jb.*(ang_t.*dellip1(:,:,3)))));
                    term1 = -ang1.*r.*st/ks.*dj_lp.*b_li;
                    term2 = r.^3/ks*lp*(lp+1).*st.*j_lp.*p_lpmp.*b_li.*(ellip1.*tau_limi-ellip2.*pi_limi*1i*mi);
                    J21(ni,np) = sum(sum(weight_map.*prefactor.*(term1+term2)));
                    dJ21_r = -ang1/ks.*st.*(dj_lp.*b_li+r.*(ks*dj2_lp.*b_li+k*d1b_li.*dj_lp));
                    ang_t = lp*(lp+1)*p_lpmp.*tau_limi/ks.*st;
                    ang_p = -1i*lp*(lp+1)*mi*p_lpmp.*pi_limi/ks.*st;
                    dJ21(ni,np,1) = sum(sum(weight_map.*prefactor.*((dJ21_r+(ang_t.*ellip1+ang_p.*ellip2).*d_r3jb).*radius_grad(:,:,1)+r3jb.*(ang_t.*dellip1(:,:,1)+ang_p.*dellip2(:,:,1)))));
                    dJ21(ni,np,2) = sum(sum(weight_map.*prefactor.*((dJ21_r+(ang_t.*ellip1+ang_p.*ellip2).*d_r3jb).*radius_grad(:,:,2)+r3jb.*(ang_t.*dellip1(:,:,2)+ang_p.*dellip2(:,:,2)))));
                    dJ21(ni,np,3) = sum(sum(weight_map.*prefactor.*((dJ21_r+(ang_t.*ellip1+ang_p.*ellip2).*d_r3jb).*radius_grad(:,:,3)+r3jb.*(ang_t.*dellip1(:,:,3)))));
                end
            end
        end
    end
end
end
        

function selection_rules = selection_rules(l,m,lp,mp,diag_switch)
if diag_switch == 1
    selection_rules = (-1)^(m)*(1+(-1)^(mp-m))*(1+(-1)^(lp+l+1));
elseif diag_switch == 2
    selection_rules = (-1)^(m)*(1+(-1)^(mp-m))*(1+(-1)^(lp+l));
else
    disp('not supported, only valid inputs for diag_switch are 1,2');
end
end
                    
                
                