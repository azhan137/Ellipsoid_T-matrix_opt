%compute weights and abscissae up to N for interval a,b

function [x,w] = generate_gauss_weights_abscissae(N,a,b)
%reform bounds
N = N-1;
N1 = N+1;
N2 = N+2;
%generate intial guess
x = linspace(-1,1,N1)';
y = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*x*N/N2);

%legendre gauss vandermond matrix
L=zeros(N1,N2);
%derivative matrix
Lp=zeros(N1,N2);
%compute zeros using recursion relation and raphson method
y0=2;

while max(abs(y-y0)) > eps
    L(:,1)=1;
    Lp(:,1)=0;
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end
    
    Lp=(N2)*(L(:,N1)-y.*L(:,N2))./(1-y.^2);
    y0=y;
    y=y0-L(:,N2)./Lp;
end

w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
x = (b-a)*y/2+(b+a)/2;