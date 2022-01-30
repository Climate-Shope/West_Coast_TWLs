function [Cg]=Cgroup(T,h)
% [Cg]=Cgroup(T,h)
% computes group velocity Cg (m/s)
% 
% T:  period (s)
% h:  depth  (m)


L=LDIS(T,h);
k=(2*pi)./L;
kh2=2.*k.*h;
n=0.5*(1+(kh2./sinh(kh2)));
C=L./T;
Cg=n.*C;