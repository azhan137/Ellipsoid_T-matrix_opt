%function that returns T and dT(a,b,c,alpha) matrices in lab frame given T
%and dT (a,b,c) in the particle frame for ellipsoids, and rotation ange alpha.
%only supports rotation about z-axis (euler angle alpha).

% lmax: maximal expansion order
% T_particle: T-matrix in particle frame
% dT_particle: derivative T-matrices foraxes a,b,c for ellipsoid
% alpha: rotation to be rotated in x-y plane, corresponds to euler angle
% alpha


function [T_lab,dT_lab] = axial_rotation(lmax,T_particle,dT_particle,alpha)

T_lab = zeros(size(T_particle));
dT_lab = zeros([size(T_particle),4]);
dT_alpha = zeros(size(T_particle));


%if the angle rotated is less than 2e-3 then consider the particle not
%rotated in which case T_lab = T_particle, dT_lab = dT_particle, and we
%only need to calculate dT_alpha
if abs(alpha) < 1e-6
    T_lab = T_particle;
    dT_lab(:,:,1:3) = dT_particle;
    
    %compute angular derivative of T-matrix
    for ti = 1:2
        for tp = 1:2
            for li = 1:lmax
                for lp = 1:lmax
                    for mi = -li:li
                        ni = multi2single_index(1,ti,li,mi,lmax);
                        for mp = -lp:lp
                            np = multi2single_index(1,tp,lp,mp,lmax);
                            dT_alpha(ni,np) = -1i*(mi-mp)*exp(-1i*alpha*(mi-mp))*T_particle(ni,np);
                        end
                    end
                end
            end
        end
    end
    dT_lab(:,:,4) = dT_alpha;
%if the particle rotation angle is nonzero, then we need to change the
%T-matrix and dT matrices
else
    %compute angular derivative of T-matrix, and rotate T-matrix and the dT-matrices 
    %from the particle frame to lab frame
    for ti = 1:2
        for tp = 1:2
            for li = 1:lmax
                for lp = 1:lmax
                    for mi = -li:li
                        ni = multi2single_index(1,ti,li,mi,lmax);
                        for mp = -lp:lp
                            np = multi2single_index(1,tp,lp,mp,lmax);
                            T_lab(ni,np) = exp(-1i*alpha*(mi-mp))*T_particle(ni,np);
                            dT_lab(ni,np,1:3) = exp(-1i*alpha*(mi-mp))*dT_particle(ni,np,:);
                            dT_lab(ni,np,4) = -1i*(mi-mp)*exp(-1i*alpha*(mi-mp))*T_particle(ni,np);
                        end
                    end
                end
            end
        end
    end
end