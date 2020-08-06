clc;
clear all;
close all;

lmax = 4; %maximum coordinate
a = 400; %along x
b = 300; %along y
c = 500; %along z
alpha = pi/3; %orientation angle, alpha as conventional in Euler angle rotations
Ntheta = 40; %number of samples for surface along polar angle
Nphi = 40; %number of samples for surface along radial angle
ni = 1; %index of medium
ns = 3; %index of scatterer
lambda = 1000;
%step sizes of radii in nm to check derivatives wrt to axes a (x), b (y), c (z)
step_size = [1e1,1,1e-1,1e-2,1e-3,1e-4];
%step sizes of angles in radians to check derivative wrt to angle
rad_size = pi./[1e1,1e2,1e3,1e4,1e5];

tic
%compute the T matrix (T_0), and derivative matrices dT_0 (1, 2, 3, 4)
%1 dT/da, 2 dT/db, 3 dT/dc, 4 dT/d(alpha)
[T_0, dT_0] = compute_T(lmax,Ntheta,Nphi,[a,b,c,alpha],ni,ns,lambda);
%perform rotation to get correct orientation
[T_0, dT_0] = axial_rotation(lmax,T_0,dT_0,alpha);
toc

%see the T-matrix
figure
imagesc(abs(T_0))
colorbar

%see the derivatives
figure
subplot(2,2,1)
imagesc(abs(squeeze(dT_0(:,:,1))))
axis equal tight
title('a')
colorbar
subplot(2,2,2)
imagesc(abs(squeeze(dT_0(:,:,2))))
axis equal tight
title('b')
colorbar
subplot(2,2,3)
imagesc(abs(squeeze(dT_0(:,:,3))))
axis equal tight
title('c')
colorbar
subplot(2,2,4)
imagesc(abs(squeeze(dT_0(:,:,4))))
axis equal tight
title('alpha')
colorbar


% %% compute error along axis a. Analytical vs Numerical
% for i = 1:length(step_size)
%     T_i = compute_T(lmax,Ntheta,Nphi,[a+step_size(i),b,c,alpha],ni,ns,lambda);
%     numerical_derivative = abs(T_i-T_0)/step_size(i);
% %     figure
% %     subplot(2,1,1)
% %     imagesc(numerical_derivative)
% %     title(strcat('numerical derivative ',num2str(step_size(i))))
% %     colorbar
% %     subplot(2,1,2)
% %     imagesc(abs(numerical_derivative-abs(squeeze(dT_0(:,:,1)))))
% %     title('error')
% %     colorbar
%     
%     max_error_a(i) = max(max(abs(numerical_derivative-abs(squeeze(dT_0(:,:,1))))));
% end
% figure
% plot(step_size, max_error_a)
% title('a error')
% xlabel('step size')
% ylabel('max error')

% %% compute error along axis b. Analytical vs Numerical
% for i = 1:length(step_size)
%     T_i = compute_T(lmax,Ntheta,Nphi,[a,b+step_size(i),c,alpha],ni,ns,lambda);
%     numerical_derivative = abs(T_i-T_0)/step_size(i);
% %     figure
% %     subplot(2,1,1)
% %     imagesc(numerical_derivative)
% %     title(strcat('numerical derivative ',num2str(step_size(i))))
% %     colorbar
% %     subplot(2,1,2)
% %     imagesc(abs(numerical_derivative-abs(squeeze(dT_0(:,:,2)))))
% %     title('error')
% %     colorbar
%     
%     max_error_b(i) = max(max(abs(numerical_derivative-abs(squeeze(dT_0(:,:,2))))));
% end
% figure
% plot(step_size, max_error_b)
% title('b error')
% xlabel('step size')
% ylabel('max error')


% %% compute error along axis c. Analytical vs Numerical
% for i = 1:length(step_size)
%     T_i = compute_T(lmax,Ntheta,Nphi,[a,b,c+step_size(i),alpha],ni,ns,lambda);
%     numerical_derivative = abs(T_i-T_0)/step_size(i);
% %     figure
% %     subplot(2,1,1)
% %     imagesc(numerical_derivative)
% %     title(strcat('numerical derivative ',num2str(step_size(i))))
% %     colorbar
% %     subplot(2,1,2)
% %     imagesc(abs(numerical_derivative-abs(squeeze(dT_0(:,:,3)))))
% %     title('error')
% %     colorbar
%     
%     max_error_c(i) = max(max(abs(numerical_derivative-abs(squeeze(dT_0(:,:,3))))));
% end

% figure
% plot(step_size, max_error_c)
% title('c error')
% xlabel('step size')
% ylabel('max error')


%% compute error along angle alpha. Analytical vs Numerical
% for i = 1:length(rad_size)
%     T_i = compute_T(lmax,Ntheta,Nphi,[a,b,c,alpha+rad_size(i)],ni,ns,lambda);
%     [T_i ~] = axial_rotation(lmax,T_i,dT_0(:,:,1:3),alpha+rad_size(i));
%     numerical_derivative = abs(T_i-T_0)/rad_size(i);
% %     figure
% %     subplot(2,1,1)
% %     imagesc(numerical_derivative)
% %     title(strcat('numerical derivative ',num2str(rad_size(i))))
% %     colorbar
% %     subplot(2,1,2)
% %     imagesc(abs(numerical_derivative-abs(squeeze(dT_0(:,:,4)))))
% %     title('error')
% %     colorbar
%     
%     max_error_rad(i) = max(max(abs(numerical_derivative-abs(squeeze(dT_0(:,:,4))))));
% end

% figure
% plot(rad_size, max_error_rad)
% title('angular error')
% xlabel('step size')
% ylabel('max error')
