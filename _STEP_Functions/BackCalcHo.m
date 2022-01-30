function [Ho,HoTwo]=BackCalcHo(T,h,H)
% function [Ho]=BackCalcHo(T,h,H)
% function to back- calcualte Ho with H at given depth 
% T: period (s)
% h: water depth (m)
% H: wave height at h (m)
% based on conservation of energy where (ECg)at DW=(ECg)at any depth
% resulting eqn is: Ho=Hd*sqrt((2*Cgd)/1.56T) where d denotes 'at given depth'

for j=1:length(h)
    Cg=Cgroup(T(j,1),h(j,1));
    Co=1.56*T(j,1);
    [Ks,kh]=shoal(T(j,1),h(j,1));  % shoaling coefficient
    Ho(j,1)=H(j,1)/Ks;
end



% second method
g=9.81;
for j=1:length(h)
    L=LDIS(T(j,1),h(j,1));
    k=(2*pi)/L;
    kh2=2*k*h(j,1);
    n=0.5*(1+(kh2/sinh(kh2)));
    C=L/T(j,1);
    HoTwo(j,1)=H(j,1)*sqrt(C*n*((4*pi)/(g*T(j,1))));
end

