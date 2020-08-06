%returns chain term of differential area element
%index 1-3, terms for a-c

function [dellip1,dellip2] = ellipsoid_function_deriv(st,ct,sp,cp,a,b,c)

dellip1(:,:,1) = -2.*cp.^2/a^3.*st.*ct;
dellip1(:,:,2) = -2.*sp.^2/b^3.*st.*ct;
dellip1(:,:,3) = 2/c^3.*st.*ct;
dellip1(:,:,4) = 2*sp.*cp.*(1/b^2-1/a^2).*st.*ct;
%dellip1(:,:,4) = (ct.^3.*st.^3.*(cp.^2/a^2+sp.^2/b^2)* ...
%    (-1/c^2+cp.^2/a^2+sp.^2/b^2))/c^2;

dellip2(:,:,1) = 2/a^3.*sp.*cp.*st;
dellip2(:,:,2) = -2/b^3.*sp.*cp.*st;
dellip2(:,:,4) = (-1/a^2+1/b^2).*st.*(cp.^2-sp.^2);
%dellip2(:,:,4) = ((-1/a^2+1/b^2)*ct.^2*cp.*st.^3*sp*(cp.^2/a^2+sp.^2/b^2))/c^2;