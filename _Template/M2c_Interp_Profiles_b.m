warning off all 

%%Second step for elevation profile interpolation

%% CONFIG 
NAME_GRD_transects  = 'Douglas'; 
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

%% 
VISIBILITY = 'on'

%% directories - data and results 
dirout      = ['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_profiles_wave_model']; if ~exist([dirout],'dir'), mkdir(dirout); end

%% profile parameters 
h1 = 0; 
h2 = +55; % profile limits, bathymetry <0 
resolution_m = 1; % m 

prof_cut_len=50;% The length that a series of continuous points must be above cutoff for it to not be discarded as a small island or rock

%% READ DATA 

file_transects          = ['E:\West_Coast_TWL_Hazards\01_Data\Transects' filesep NAME_GRD_transects '_Transects_v3.mat']; 
file_dem_utm            = ['E:\West_Coast_TWL_Hazards\01_Data\DEM' filesep NAME_GRD_transects '_Coastal_DEM.tif']; 
%MHW level and MSL level by transect
file_MHW                =['E:\West_Coast_TWL_Hazards\03_Results\MHW' filesep NAME_GRD_transects '_MHW.mat']; 
file_MSL                =['E:\West_Coast_TWL_Hazards\03_Results\MSL' filesep NAME_GRD_transects '_MSL.mat']; 

%%
num_dem=1; %User defined beased on avaialbe DEMs
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
Ntr             = numel(transects.X(:,2)); 

% load R-Topo and Fbound meta information 
for rr=1:num_dem;
    load(['Rtopo_' num2str(rr)]);
    load(['F_bound_' num2str(rr)]);
end

[PATHSTR2,NAME2,EXT]= fileparts(file_dem_utm);

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

    %need to make the square from the worldlimits
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




% 2) FIND CORRECT GRID 
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



depth_lim = [h1 h2]; % profile limits, bathymetry <0 

array_for_plots = [0]; 

time1 = now; 

transects_mod = transects; 
transects_mod.xoff=transects.X(:,1);
transects_mod.yoff=transects.Y(:,1);
transects_mod.xon=transects.X(:,3);
transects_mod.yon=transects.Y(:,3);
transects_mod.x0=transects.X(:,2);
transects_mod.y0=transects.Y(:,2);
INCORRECT_RUNS = [];


for ii = 1:Ntr;
    OBJECTID = transects.OBJECTID(ii); 
    disp([num2str(ii),'/',num2str(Ntr)]) 
    

    Lx = transects_mod.xoff(ii)-transects_mod.xon(ii); 
    Ly = transects_mod.yoff(ii)-transects_mod.yon(ii); 
    
    L0 = sqrt(Lx.^2 + Ly.^2);  % true distance, from UTMs 
        
    Npts = L0./resolution_m; 
    Npts = round(Npts); 
    
    %get the first intersection from X0,Y0 to the boundary
        
    xpf   = linspace(transects_mod.xoff(ii) ,transects_mod.xon(ii) ,Npts);
    ypf   = linspace(transects_mod.yoff(ii) ,transects_mod.yon(ii) ,Npts);
    
    [~,xt] = min(abs(xpf-transects_mod.x0(ii)));
    [~,yt] = min(abs(ypf-transects_mod.y0(ii)));
    

    uu=GRIDid(ii);
    eval(['temp=F_bound_' num2str(uu) ';']);
   
    L1=[xpf;ypf];
    if length(temp.X(:,1))>1;
       L2=[temp.X';temp.Y'];
    else
    L2=[temp.X;temp.Y];
    end
    P = InterX(L1,L2);
   
    if length(P(1,:))>1;
    xi=[];
    yi=[];
    for yy=1:length(P(1,:));
    [~,xP] = min(abs(xpf-P(1,yy)));
    [~,yP] = min(abs(ypf-P(2,yy)));
   
    xi=[xi, xP];
    yi=[yi, yP];
    end
   
    xi=sort(xi);
    yi=sort(yi);
    
    %find closest to xt, but greater 
    idx=xi;
    idx(idx<xt)=[];
    [~,idx2]=min(abs(idx-xt));
    
    xend=idx(idx2);
    yend=idx(idx2);
    
    xpf(idx+1:end)=[];
    ypf(idx+1:end)=[];
    
    %find the closest point in the series that is not greater
    xi(xi>xt)=[];
    yi(yi>yt)=[];
    
    %end in series should be the edge
    xi=xi(end);
    yi=yi(end);
    
    xpf(1:xi-1)=[];
    ypf(1:xi-1)=[];
    else
    
     in=inpolygon(xpf,ypf,temp.X,temp.Y);
     in=find(in);
     xpf=xpf(in);
     ypf=ypf(in);
    end
    
    
    if ~isempty(xpf);
    Transects2.xpf(ii,:)=[xpf(1) xpf(end)];
    Transects2.ypf(ii,:)=[ypf(1) ypf(end)];
    else
    Transects2.xpf(ii,:)=[NaN NaN];
    Transects2.ypf(ii,:)=[NaN NaN]; 
    end
    
end

save('Transects2','Transects2');