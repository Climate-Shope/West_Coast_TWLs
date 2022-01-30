function Hb=Hbreaking(Ho,T,gammab)
% function to calculate the breaking wave height (from conservation of energy)
% b=Hbreaking(Ho,T,gammab)

g=9.81;
Hb=(Ho.^4.*((g.*T)/(4.*pi)).^2.*(gammab/g)).^(1/5);