%Interpolate profile elevation from DEM at cross shore transects

warning off all 
%% CONFIG 

NAME_GRD_transects  = 'Douglas'; 
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
%% 
VISIBILITY = 'on'

%% 
dirout      = ['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_profiles_wave_model']; if ~exist([dirout],'dir'), mkdir(dirout); end

%% profile parameters 
h1 = 0; 
h2 = +55; % profile limits, bathymetry <0 
resolution_m = 1; % m 

%profile sizing parameters
prof_cut_len=50;% The length that a series of continuous points must be above cutoff for it to not be discarded as a small island or rock

%% READ DATA 

file_transects          = ['E:\West_Coast_TWL_Hazards\01_Data\Transects' filesep NAME_GRD_transects '_Transects_v3.mat']; 
file_dem_utm            = ['E:\West_Coast_TWL_Hazards\01_Data\DEM' filesep NAME_GRD_transects '_Coastal_DEM.tif']; 
file_MHW                =['E:\West_Coast_TWL_Hazards\03_Results\MHW' filesep NAME_GRD_transects '_MHW.mat']; 
file_MSL                =['E:\West_Coast_TWL_Hazards\03_Results\MSL' filesep NAME_GRD_transects '_MSL.mat']; 


%% number of DEMs used for intepolation
num_dem=1;
%%
load(file_transects);
eval(['transects =' NAME_GRD_transects '_Transects_v3;']);

load(file_MHW);
load(file_MSL);

if isfield(transects, 'OBJECTID')==0
        for ii = 1:numel(transects.X(:,2)) 
            transects.OBJECTID(ii) = ii; 
        end
end
Ntr  = numel(transects.X(:,2)); 

for rr=1:num_dem;
    
    if num_dem==1;
       file_dem_utm         = ['E:\West_Coast_TWL_Hazards\01_Data\DEM' filesep NAME_GRD_transects '_Coastal_DEM_1m.tif'];
    else
        
        file_dem_utm            = ['E:\West_Coast_TWL_Hazards\01_Data\DEM' filesep NAME_GRD_transects '_' num2str(rr) '_Coastal_DEM_1m.tif']; 
    end


%Get the meta infomation for the DEM
eval(['[topo, Rtopo_' num2str(rr) '] = geotiffread(file_dem_utm);']);
clear topo
end


%Save the Rtopo Data
for rr=1:num_dem;
    save(['E:\West_Coast_TWL_Hazards\_STEP\_' NAME_GRD_transects filesep 'Rtopo_' num2str(rr)],['Rtopo_' num2str(rr)]);
end



[PATHSTR2,NAME2,EXT]= fileparts(file_dem_utm);



%--------------------------------------------------------------------------
% CHECK PROJECTION OF EACH RASTER!!!!
%figure, hold on 
%plot([Rtopo.XLimWorld(:); flipud(Rtopo.XLimWorld(:))],[Rtopo.YLimWorld(1); Rtopo.YLimWorld(1); Rtopo.YLimWorld(2); Rtopo.YLimWorld(2)],'-k')

%--------------------------------------------------------------------------

%% FIND SWAN GRID FOR EACH TRANSECT 
SWANgrdcenter=zeros(num_dem,2); 

for grid =1:num_dem;
    
    eval(['zopo=Rtopo_' num2str(grid)]);
    X=zopo.XWorldLimits;
    Y=zopo.YWorldLimits;
    SWANgrdcenter(grid,1)=mean([min(X(:)) max(X(:))]); 
    SWANgrdcenter(grid,2)=mean([min(Y(:)) max(Y(:))]); 
    
    
    testx=zopo.XWorldLimits(1)+1:zopo.XWorldLimits(end);
    testy=zopo.YWorldLimits(1):zopo.YWorldLimits(end);
    testy=flipud(testy');
    
end

%% find which grid to use 
disp('finding right grids...') 

coords = [transects.X(:,2) transects.Y(:,2)];

Ntr     = numel(coords(:,1));  

DIST =[]; 

clear xv yv
INPOINTS = nan(numel(coords(:,1)),num_dem); 

for goose = 1:num_dem;
    eval(['zopo=Rtopo_' num2str(goose)]);
    X=zopo.XWorldLimits;
    Y=zopo.YWorldLimits;
    
    xx=[X(1); X(2); X(2); X(1)];
    yy=[Y(1); Y(1); Y(2); Y(2)];
    
    figure, plot(xx,yy, '.-')
    hold on, 
    plot(coords(:,1),coords(:,2),'.') 
    in = inpolygon(coords(:,1),coords(:,2),xx,yy);    
    INPOINTS(:,goose)=in; 
    plot(coords(in,1),coords(in,2),'.r')
end

GRIDid = nan(Ntr,1); 
disp('finding right grid for each point...') 
temp=INPOINTS;
for ii =1:num_dem;
    temp(:,ii)=temp(:,ii)*ii;
end

for ii=1:num_dem;
    cal=find(INPOINTS(:,ii));
    GRIDid(cal)=ii;
end

GridL=unique(GRIDid);
temp=INPOINTS;
for hh=1:num_dem-1;
    bb=GridL(hh);
    bb=find(temp(:,hh)>0);
    startpoint=bb(end);
    bb=GridL(hh+1);
    bb=find(temp(:,hh+1)>0);
    endpoint=bb(1);
    d=round((startpoint-endpoint)/2);
    e=(endpoint):(endpoint+d);
    GRIDid(e)=GridL(hh);
    e=(startpoint-d+1):(startpoint);
    GRIDid(e)=GridL(hh)+1;
    

end








figure, hold on 
plot(coords(:,1), coords(:,2),'.k') 
plot(coords(GRIDid==1,1), coords(GRIDid==1,2),'.r') 
plot(coords(GRIDid==2,1), coords(GRIDid==2,2),'.g') 
plot(coords(GRIDid==3,1), coords(GRIDid==3,2),'.m') 
plot(coords(GRIDid==4,1), coords(GRIDid==4,2),'.y') 
plot(coords(GRIDid==5,1), coords(GRIDid==5,2),'.b') 
plot(coords(GRIDid==6,1), coords(GRIDid==6,2),'.','color',[0.2 0.5 0.9]) 

temp2=find(isnan(GRIDid));
temp3=find(~isnan(GRIDid));
for tt=1:length(temp2);
    pos=temp2(tt);
    newgrid=pos-1;
    if newgrid==0;
        newgrid=temp3(1);
    end
    GRIDid(temp2(tt))=GRIDid(newgrid);
end


GridL=unique(GRIDid);
%Now calculate Scattered interpolant variable F, which is saved separately,
%for the next set of interpolation scripts
for kk=1:num_dem;
    eval(['zopo=Rtopo_' num2str(kk) ';']);



    if num_dem==1;
       file_dem_utm            = ['E:\West_Coast_TWL_Hazards\01_Data\DEM' filesep NAME_GRD_transects '_Coastal_DEM_1m.tif'];
    else
        
    file_dem_utm            = ['E:\West_Coast_TWL_Hazards\01_Data\DEM' filesep NAME_GRD_transects '_' num2str(kk) '_Coastal_DEM_1m.tif']; 
    end

[topo, ~] = geotiffread(file_dem_utm);
topo(topo<-1e10) = NaN; 


A=~isnan(topo);

[row col]=find(A);


Z=[];
for ii=1:length(row);
    try
   
    Z(ii)=double(topo(row(ii), col(ii)));
    catch
     Z(ii)=NaN;
    end
end
clear topo




[XX,YY]= meshgrid(linspace(zopo.XLimWorld(1), zopo.XLimWorld(2),zopo.RasterSize(2)), linspace(zopo.YLimWorld(1), zopo.YLimWorld(2),zopo.RasterSize(1))); 
YY = flipud(YY); 

X=[];
Y=[];
for ii=1:length(row);
    try
   
    X(ii)=XX(row(ii), col(ii));
    Y(ii)=YY(row(ii), col(ii));
    catch
     X(ii)=NaN;
     Y(ii)=NaN;
    end
end
clear XX YY

temp1=find(isnan(Z));
X(temp1)=NaN;
Y(temp1)=NaN;

temp1=find(isnan(X));
Z(temp1)=NaN;

Z(isnan(Z))=[];
X(isnan(X))=[];
Y(isnan(Y))=[];



F = scatteredInterpolant(X',Y',Z','nearest','none');
temp=boundary(X',Y');
F_bound.X=X(temp);F_bound.Y=Y(temp);%Bounding polygon of the DEM
clear X Y Z 

eval(['F_' num2str(kk) '=F;']);
eval(['F_bound_' num2str(kk) '=F_bound;']);
save(['E:\West_Coast_TWL_Hazards\_STEP\_' NAME_GRD_transects filesep 'F_' num2str(kk)],['F_' num2str(kk)]);
save(['E:\West_Coast_TWL_Hazards\_STEP\_' NAME_GRD_transects filesep 'F_bound_' num2str(kk)],['F_bound_' num2str(kk)]);
eval(['clear F_' num2str(kk)]);
eval(['clear F_bound_' num2str(kk)]);

end































