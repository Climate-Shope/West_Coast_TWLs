function [Ks,kh]=shoal(T,h)
% [Ks,kh]=shoal(T,h)

g=9.81;


% Deep water wave characterstics
Co=(g*T)/(2*pi);

% local water depth wave characterstics
L=LDIS(T,h);
C=L/T;
k=(2*pi)/L;
kh2=2*k*h;
Cg=C*(1+kh2/sinh(kh2))/2;
Ks=sqrt(Co/(2*Cg));
kh=k*h;
