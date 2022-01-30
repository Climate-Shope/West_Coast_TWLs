function [xi,yi] = Get_Spaced_Points(Xutm, Yutm, delta);
%[xi,yi] = Get_Spaced_Points(Xutm, Yutm, delta, buffer);
% The function Get_Spaced_Points marchs along a shorleine and every
% delta (in meters) it calculates the location of the shorleine at that
% point. The shoreline should be relatively edited first and the vectors of points in
% order (preferably marching from south to north for future calculations)

%Discontinuous data is okay as the function will aply points over linearly
%interpolated lines within gaps. 

%%% If there is a domain polygon, this will only keep the points within the
%%% domain, so do not have to worry about estuaries and harbors in data
%%% gaps.

%delta=10;%meters
X_comp=(diff(Xutm)).^2;
Y_comp=(diff(Yutm)).^2;
Dist=sqrt(X_comp+Y_comp);
Cum_Dist=sum(Dist);
npts=Cum_Dist/delta;
npts=floor(npts);

%% Now have the number of points needed at more or less 10m resolution
%% Interpolate Xutm and Yutm to those points

Cum_Dist=cumsum(Dist);
Cum_Dist=[0;Cum_Dist];

dx=Cum_Dist(end)/npts;
xi=interp1(Cum_Dist,Xutm,Cum_Dist(1):delta:Cum_Dist(end),'linear');
yi=interp1(Cum_Dist,Yutm,Cum_Dist(1):delta:Cum_Dist(end),'linear');

%There will be variations in the distances because it is generated
%every 10m alongshore, which when meauring the coarser distances, loses some
%of the detail, so points may appear closer due to shoreline crenulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


