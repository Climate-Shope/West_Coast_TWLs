addpath(['E:\West_Coast_TWL_Hazards\_STEP' filesep '_SRC_all_runs\RBF_EscalaresDireccionales']) 
set(0, 'DefaultFigureColor','White','DefaultFigurePaperPositionMode','auto')
%%  reconstruct time series from the Swan progations Part A to set up for use with another machine
%   Adapted from the WAVE DATA FOR REEFS - WD4R project
% 
%% CONFIG 
NAME_GRD_gow        = {'NAWC29','NAWC29'};
NAME_GRD_transects  = 'Douglas';  
NAME_GRD_swan       = {'wavm-hw12-LDUTM_100m_1','wavm-hw12-LDUTM_100m_2'};
NAME_GRD_Bathy       ={'LDUTM_100m_1','LDUTM_100m_2'};
DirBathy='E:\West_Coast_TWL_Hazards\01_Data\SWAN_Grids\Oregon_All_grids';

%% First Determine when ot Use each GOW output point 

load(['E:\West_Coast_TWL_Hazards\03_Results\Transects\Offshore_15m_Locations' filesep NAME_GRD_transects filesep NAME_GRD_transects '_Contour_pts.mat']);
%% FIND SWAN GRID FOR EACH TRANSECT 
SWANgrdcenter=zeros(numel(NAME_GRD_swan),2); 
disp('reading swan grids...') 
clear SWAN GRID SWANgrd*

for grid =1:numel(NAME_GRD_Bathy)
    
    FILENAME = [DirBathy filesep NAME_GRD_Bathy{grid} '.grd'];
    [X,Y,ENC] = wlgrid('read', FILENAME);


    SWANgrdcenter(grid,1)=mean([min(X(:)) max(X(:))]); 
    SWANgrdcenter(grid,2)=mean([min(Y(:)) max(Y(:))]); 
    
    SWAN{grid}.X= X;    
    SWAN{grid}.Y= Y;     
end
%% find which grid to use 
disp('finding right grids...') 


eval(['head_transects=' NAME_GRD_transects '_Contour_pts;']);
coords = [head_transects.X' head_transects.Y'];

Ntr     = numel(coords(:,1)); 

DIST =[]; 

clear xv yv
INPOINTS = nan(numel(coords(:,1)),numel(NAME_GRD_Bathy)); 
for grid = 1:numel(NAME_GRD_swan)
    X = SWAN{grid}.X; 
    Y = SWAN{grid}.Y; 
    
   
    
     points=[X(:),Y(:)];
    
    points(isnan(points(:,1)),:)=[];
    
    k=boundary(points(:,1),points(:,2));
    
    xx=points(k,1);
    yy=points(k,2);

    
    figure, plot(xx,yy, '.-')
    hold on, 
    plot(coords(:,1),coords(:,2),'.') 
    in = inpolygon(coords(:,1),coords(:,2),xx,yy);    
    INPOINTS(:,grid)=in; 
    plot(coords(in,1),coords(in,2),'.r')
end

% 2) FIND RIGHT GRID 
GRIDid = nan(Ntr,1); 
disp('finding right grid for each point...') 
for ii =1:Ntr
    VALUE =nan(numel(NAME_GRD_swan),1); row=VALUE; col=VALUE; 
    for grid = 1:numel(NAME_GRD_swan)
        if INPOINTS(ii,grid)==0, 
            continue, 
        else
            [VALUE(grid), row(grid),col(grid)] = nearestpoint_coords(SWAN{grid}.X,SWAN{grid}.Y, coords(ii,1),coords(ii,2)); 
        end % outside the grid
        
    end
    if 1-all(isnan(VALUE))
        [temp,GRIDid(ii)]=nanmax(VALUE); % point with largest distance to contour of grid 
    else 
       GRIDid(ii)=nan; 
    end
end


GridL=unique(GRIDid);



%%
dirout =['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\time_series_reconstructed'];
if ~exist([pwd filesep dirout],'dir'), mkdir(dirout); end

%% LOAD TRANSECTS 

eval(['Ntr=length(' NAME_GRD_transects '_Contour_pts.X);']);
%% PROPAGATE WAVES - PARAMETERS 

% PARAMETERS FOR PROPAGATION 
load(['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_swan\' NAME_GRD_transects '_Hsinterp.mat']);
load(['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_swan\' NAME_GRD_transects '_hwinterp.mat']);
load(['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_swan\' NAME_GRD_transects '_Tminterp.mat']);
load(['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_swan\' NAME_GRD_transects '_Tpinterp.mat']);



%% load classification for grid 
SIZE = numel(Hsinterp(:,1)); %500; 
rad=pi/180; 

for qq=1:numel(NAME_GRD_gow);
% DIRECTORY WITH THE CLASSIFICATION - FROM M1 CODES 

directories.classification = ['E:\West_Coast_TWL_Hazards\03_Results\results_GOW_3',filesep,NAME_GRD_gow{qq}]; 
%500 wave regime classifications form GOW data provided from Reguero et al.
%2012
classification  = load([directories.classification filesep 'MDA_Bmus_' NAME_GRD_gow{qq} '_' num2str(SIZE) '.mat']) 
positions       = load([directories.classification filesep 'time_and_positions_' NAME_GRD_gow{qq} '_' num2str(SIZE) '_MDA' '.mat']);

%% load X matrix for reconstruction with RBF 
%GOW data provided from Reguero et al. 2012
ALLDATA = load([directories.classification filesep 'X_' NAME_GRD_gow{qq},'.mat']);

%% CHECK UNIQUE COLUMNS IN CLASSIFICATION 
IND = find(unique(classification.final(1,:))); 
if numel(IND)~=numel(classification.final(1,:)), 
%     error('no unique columns') 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    REMCOLS = 10:12; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    classification.subset(:,REMCOLS)=[]; 
    classification.final(:,REMCOLS)=[]; 
    classification.escalar(:,[7,8])=[]; 
    classification.direccional(:,4)=[]; 
    
    
    ALLDATA.X(:,REMCOLS)=[]; 
    ALLDATA.datos_n(:,REMCOLS)=[]; 
    ALLDATA.labels(REMCOLS)=[]; 
    ALLDATA.escalar = classification.escalar; 
    ALLDATA.direccional = classification.direccional; 
    
end

%% Okay, save each output as a separate piece
eval(['ALLDATA_' num2str(qq) '=ALLDATA;']);
eval(['classification_' num2str(qq) '=classification;']);
eval(['positions_' num2str(qq) '=positions;']);
end







%%
load(['E:\West_Coast_TWL_Hazards\01_Data\GOW_Data',filesep,'time_series.mat'])  % loads time 
timevec = datevec(time); 
indsremove = (timevec(:,1)==1948 & timevec(:,2)==1); 
timevec(indsremove,:)=[]; 
time = datenum(timevec); 

fname= [dirout,filesep,'time']; 
m=matfile(sprintf('%s.mat', fname),'writable',true); 
m.time = time; 
m.timevec  = timevec; 
clear fname 

%% PROPAGATE WAVES - PARAMETERS 

% PARAMETERS FOR PROPAGATION 
load(['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_swan\' NAME_GRD_transects '_Hsinterp.mat']);
ii_grid_array = find(1-isnan(Hsinterp(1,:))); % only transects in SWAN grid 
array_transects = 1:Ntr;  % full array of transects. 
Nact_trs= numel(ii_grid_array); 

eval(['save ''',dirout,filesep,'ii_grid_array.mat'' ii_grid_array Nact_trs'])

% Ntr
%% for every transect
t1 = now;
MAKE_PLOT = 1; 


%% save workspace 

save('Douglas_TS_Input');
