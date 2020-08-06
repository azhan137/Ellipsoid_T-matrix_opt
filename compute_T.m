%%Written and implemented by Alan Zhan
%No guarantee the code is in any working order, and no warranty is implied.



%%compute T-matrix of an ellipsoid
%ellip_params: 3x1 matrix: a,b,c along x,y,z
%ni: refractive index of surrounding medium
%ns: refractive index of particle
%lambda: incident wavelength
%Nphi, Ntheta: number of degrees for quadrature for phi, theta respectively
%lmax: maximum multipole order

%see pc waterman paper on extended boundary condition method
%returns T and dT in particle frame

function [T,dT] = compute_T(lmax,Ntheta,Nphi,ellip_params,ni,ns,lambda)
a = ellip_params(1,1);
b = ellip_params(1,2);
c = ellip_params(1,3);
% 
[Q,dQ] = compute_Q(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,3);
[rQ,drQ] = compute_Q(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,1);
Qinv = inv(Q);
T_particle = rQ*Qinv;
dT_particle = zeros([size(T_particle),3]);
dT_particle(:,:,1) = (drQ(:,:,1)-T_particle*dQ(:,:,1))*Qinv;
dT_particle(:,:,2) = (drQ(:,:,2)-T_particle*dQ(:,:,2))*Qinv;
dT_particle(:,:,3) = (drQ(:,:,3)-T_particle*dQ(:,:,3))*Qinv;

T = T_particle;
dT = dT_particle;

% alpha = ellip_params(1,4);
% T_lab = zeros(size(T_particle));
% dT_lab = zeros([size(T_particle),4]);
% dT_alpha = zeros(size(T_particle));

% %if the angle rotated is less than 1e-6 then consider the particle not
% %rotated in which case T_lab = T_particle, dT_lab = dT_particle, and we
% %only need to calculate dT_alpha
% if abs(alpha) < 1e-6
%     T_lab = T_particle;
%     dT_lab(:,:,1:3) = dT_particle;
%     
%     %compute angular derivative of T-matrix
%     for ti = 1:2
%         for tp = 1:2
%             for li = 1:lmax
%                 for lp = 1:lmax
%                     for mi = -li:li
%                         ni = multi2single_index(1,ti,li,mi,lmax);
%                         for mp = -lp:lp
%                             np = multi2single_index(1,tp,lp,mp,lmax);
%                             dT_alpha(ni,np) = -1i*(mi-mp)*exp(-1i*alpha*(mi-mp))*T_particle(ni,np);
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     dT_lab(:,:,4) = dT_alpha;
% %if the particle rotation angle is nonzero, then we need to change the
% %T-matrix and dT matrices
% else
%     %compute angular derivative of T-matrix, and rotate T-matrix and the dT-matrices 
%     %from the particle frame to lab frame
%     for ti = 1:2
%         for tp = 1:2
%             for li = 1:lmax
%                 for lp = 1:lmax
%                     for mi = -li:li
%                         ni = multi2single_index(1,ti,li,mi,lmax);
%                         for mp = -lp:lp
%                             np = multi2single_index(1,tp,lp,mp,lmax);
%                             T_lab(ni,np) = real(exp(-1i*alpha*(mi-mp)))*T_particle(ni,np);
%                             dT_lab(ni,np,1:3) = real(exp(-1i*alpha*(mi-mp)))*dT_particle(ni,np,:);
%                             dT_lab(ni,np,4) = 1i*real(1i*(mi-mp)*exp(-1i*alpha*(mi-mp)))*T_particle(ni,np);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

% T = T_lab;
% dT = dT_lab;