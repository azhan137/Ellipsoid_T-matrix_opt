
%% function retrieves ellipsoid surface functions
% ellip_switch:
%1. cos^2(phi)/a^2+sin^2(phi)/b^2-1/c^2
%2. -1/a^2+1/b^2
% maybe derivatives later?

function ellip_funct = ellipsoid_function(a,b,c,phi,ellip_switch)

switch ellip_switch
    case 1
        if a == b && b == c
            ellip_funct = 0;
        elseif a == b
            ellip_funct = 1/a^2-1/c^2;
        else
            ellip_funct = (cos(phi)/a).^2+(sin(phi)/b).^2-1/c^2;
        end
    case 2
        ellip_funct = -1/a^2+1/b^2;
    otherwise
        disp('value of ellip_switch not implemented, valid values 1,2');
end
        
end