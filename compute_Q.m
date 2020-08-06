function [Q,dQ] = compute_Q(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,nu)

[J11,J12,J21,J22,dJ11,dJ12,dJ21,dJ22] = compute_J_ellip_celes_posm(lmax,Ntheta,Nphi,a,b,c,ni,ns,lambda,nu);

k = 2*pi*ni/lambda;
ks = 2*pi*ns/lambda;

P = -1i*k*(ks*J21+k*J12);
R = -1i*k*(ks*J11+k*J22);
S = -1i*k*(ks*J22+k*J11);
U = -1i*k*(ks*J12+k*J21);

dP = -1i*k*(ks*dJ21+k*dJ12);
dR = -1i*k*(ks*dJ11+k*dJ22);
dS = -1i*k*(ks*dJ22+k*dJ11);
dU = -1i*k*(ks*dJ12+k*dJ21);

% figure
% subplot(2,2,1)
% imagesc(abs(P))
% title(num2str(nu));
% colorbar
% subplot(2,2,2)
% imagesc(abs(R))
% title('R')
% colorbar
% subplot(2,2,3)
% imagesc(abs(S))
% title('S')
% colorbar
% subplot(2,2,4)
% imagesc(abs(U))
% colorbar
% title('U')

Q = [P,R;S,U];
dQ = [dP,dR;dS,dU];
