%% Script to intepolate Downscaled SWAN inputs to 15 m offhsore locations
addpath(genpath('E:\West_Coast_TWL_Hazards\OPENEARTHTOOLS\delft3d_matlab')) % delft 3d toolbox 

%% CONFIG 
NAME_GRD_transects  = 'Douglas'; 

%% DEFINE GRID SWANS IN THE DOMAIN 
NAME_GRD_swan       = {'wavm-hw12-LDUTM_100m_1','wavm-hw12-LDUTM_100m_2'};
NAME_GRD_Bathy       ={'LDUTM_100m_1','LDUTM_100m_2'};
County={'LaneDouglas','LaneDouglas'};
%% directories 
DirBathy='E:\West_Coast_TWL_Hazards\01_Data\SWAN_Grids\Oregon_All_grids';
dirout = ['E:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects filesep 'interp_swan']; 
if exist(dirout)==0, mkdir(dirout),end

%% LOAD Transect Endpoints 
load(['E:\West_Coast_TWL_Hazards\03_Results\Transects\Offshore_15m_Locations' filesep NAME_GRD_transects filesep NAME_GRD_transects '_Contour_pts.mat']);
%% FIND SWAN GRID FOR EACH TRANSECT 
SWANgrdcenter=zeros(numel(NAME_GRD_swan),2); 
disp('reading swan grids...') 
clear SWAN GRID SWANgrd*

SWANgrdcenter=zeros(numel(NAME_GRD_Bathy),2); 
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
for grid = 1:numel(NAME_GRD_Bathy)
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



GridL=unique(GRIDid,'stable');


%% interpolate, one grid once at a time 

ddd=1; 

dat_name = ['E:\West_Coast_TWL_Hazards\03_Results\Completed_SWAN_Runs\Oregon\', County{ddd} , filesep ,NAME_GRD_swan{ddd},'.dat'];
%vs_use and qpread ate part of the Deltares DEFLT3D Matlab Toolbox 
fid = vs_use(dat_name);
[Times] = qpread(fid,'hsig wave height','times');


Ntimes  = numel(Times); 

Hsinterp = nan(Ntimes,Ntr); 
Tminterp = nan(Ntimes,Ntr); 
Tpinterp = nan(Ntimes,Ntr);
hwinterp = nan(Ntimes,Ntr);


for grid =1:numel(NAME_GRD_swan)
       
    dat_name = ['E:\West_Coast_TWL_Hazards\03_Results\Completed_SWAN_Runs\Oregon\', County{grid} , filesep ,NAME_GRD_swan{grid},'.dat'];
    disp(dat_name) 
    
    fid = vs_use(dat_name);
    [DataFields,Dims,NVal] = qpread(fid);

    disp('reading outputs...')
    Hs      = qpread(fid,'hsig wave height','griddata',0); 
    Tp      = qpread(fid,'relative peak wave period','griddata',0);
    Tm      = qpread(fid,'mean wave period','griddata',0);
    hw      = qpread(fid,'water depth','griddata',0);



    %% interpolate in selected points
    ind = find(GRIDid==grid); 
    
    MASK= squeeze(nanmean(Hs.Val,1)); 
    MASK(MASK~=0)=1;
    MASK(MASK==0)=nan; 
    
    %%
    figure, pc = pcolor(Hs.X, Hs.Y, squeeze(Hs.Val(1,:,:)).*MASK), set(pc, 'linestyle','none') 
    hold on, plot(coords(:,1),coords(:,2),'ok') 
    
    %%
    X = Hs.X.*MASK; 
    Y = Hs.Y.*MASK; 
    
    for ii = 1:numel(ind)
        disp(ind(ii))

        [VALUE, row,col] = nearestpoint_coords(X,Y, coords(ind(ii),1),coords(ind(ii),2)); 
        % control by distance 
        if VALUE<200 % in meters 
            Hsinterp(:,ind(ii))=squeeze(Hs.Val(:,row,col)); 
            Tminterp(:,ind(ii))=squeeze(Tm.Val(:,row,col)); 
            Tpinterp(:,ind(ii))=squeeze(Tp.Val(:,row,col)); 
            hwinterp(:,ind(ii))=squeeze(hw.Val(:,row,col)); 
            
        end
    end
end
%%
ind2 = find(isnan(Hsinterp(1,:)));
if 1-isempty(ind2), 
    Hsinterp(1,ind2); 
    warning('points with missing values'), 
    disp(ind2) 
    disp([num2str(numel(ind2)),'/',num2str(Ntr)])
end 

%% save 

save([dirout filesep, NAME_GRD_transects,'_Hsinterp.mat'],['Hsinterp'])
save([dirout filesep, NAME_GRD_transects,'_Tminterp.mat'],['Tminterp'])
save([dirout filesep, NAME_GRD_transects,'_Tpinterp.mat'],['Tpinterp'])
save([dirout filesep, NAME_GRD_transects,'_hwinterp.mat'],['hwinterp'])
    
winopen([dirout ]) 


