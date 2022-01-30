clear all; run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
% code to generate approximate shoreline locations and cast major transects
% every 100m
% from Coastal DEM data 

NAME_GRD_transects  = 'Douglas'; 


run('run_init_directories.m')


Contour_elevation=2; %Should be an elevation greater than the variation seen in water returns 2m seems to work
Consecutive_Length=200; %Number of consecutive points for the reigon to be considered a constinuous shoreline and not an offshore island or rock 
delta=10; %desired spacing of minor transects
r_on=500; %(m) how far to extend the transect onshore
r_off=1000; %(m) how far to extend the transect offshore
num_points=1; %the points on either side of the transect used to calc its shoreline angle
dirout=directories.TRANSECTS;

%% Load DEM
file_dem_utm  = [directories.DEM filesep NAME_GRD_transects '_Coastal_DEM.tif']; 
[topo, Rtopo] = geotiffread(file_dem_utm); 

%USGS WC 2016 Lidar is 0.5m, so cut to 2m for space
topo_10m=topo(1:20:end,1:20:end); %If we can push to 2m, life gets a lot easier, for now it cuts the DEM to 1m res
clear topo
topo_10m(topo_10m<-1e10) = NaN; 
topo_10m=flipud(topo_10m);
%%

%%% Create a meshgid to make contouring 
X=linspace(Rtopo.XWorldLimits(1),Rtopo.XWorldLimits(2),length(topo_10m(1,:)));
Y=linspace(Rtopo.YWorldLimits(1),Rtopo.YWorldLimits(2),length(topo_10m(:,1)));

[XX, YY]=meshgrid(X,Y);
[C,H]=contour(XX,YY,topo_10m,[Contour_elevation Contour_elevation]);

test=find(C(1,:)==Contour_elevation);
test_d=diff(test);
test_r=find(test_d > Consecutive_Length);

along_idx_start=test(test_r)+1;
along_idx_end=test(test_r+1)-1;

New_lon=[];
New_lat=[];
for ii=1:length(along_idx_start);
    x=C(1,along_idx_start(ii):along_idx_end(ii));
    y=C(2,along_idx_start(ii):along_idx_end(ii));
    
    x=fliplr(x);
    y=fliplr(y);
    
    
    New_lon=[New_lon,x];
    New_lat=[New_lat,y];
end


%Get the points at every 10 m 

Points=[New_lon', New_lat'];
Points=unique(Points,'rows','stable');

[xi,yi] = Get_Spaced_Points(Points(:,1), Points(:,2), delta);
figure;plot(xi,yi)


[TransectX, TransectY] = Get_normal_transects(xi,yi,r_on,r_off,num_points);

TransectX=TransectX(1:10:end,:);
TransectY=TransectY(1:10:end,:);


Transects.X=TransectX;
Transects.Y=TransectY;

name=[NAME_GRD_transects '_Transects'];
eval([name '=Transects;']);
save([dirout filesep name '.mat'], name);



    














