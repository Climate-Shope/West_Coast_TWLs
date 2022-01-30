function [Ho]=BackCalcHo_Limits(T,h,H)
% function [Ho]=BackCalcHo(T,h,H)
% function to back- calcualte Ho with H at given depth 
% T: period (s)
% h: water depth (m)
% H: wave height at h (m)
% based on conservation of energy where (ECg)at DW=(ECg)at any depth
% resulting eqn is: Ho=Hd*sqrt((2*Cgd)/1.56T) where d denotes 'at given depth'
% Updated by James Shope 2019

Up=[linspace(0,max(T),1000)]';


if Up(1)>0;
edges=0;
edges=[edges;(edges(end)+Up(1))/2];
for ii=2:length(Up)-1;
    temp=(Up(ii)+Up(ii+1))/2;
    edges=[edges;temp];
end
edges=[edges;Up(end)];
else
   edges=Up(1)-0.01; 
   edges=[edges;(edges(end)+Up(1))/2];
   for ii=2:length(Up)-1;
    temp=(Up(ii)+Up(ii+1))/2;
    edges=[edges;temp];
   end
edges=[edges;Up(end)];
end

for j=1:length(Up)
    [Ks(j),kh]=shoal(Up(j,1),h);  % shoaling coefficient

end

tmp = discretize(T,edges);
Ks=Ks(tmp);


Ho=H./Ks';



