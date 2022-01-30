% Calculate offshore points to find location of the 15m contour along shore normal transects to retrive
% downscaled wave data
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
clear all
warning off all 
%% CONFIG 
NAME_GRD_transects  = 'Douglas'; 
 
VISIBILITY = 'on'

%% directories - data and results 

dirout      = ['E:\West_Coast_TWL_Hazards\03_Results\Transects\Offshore_15m_Locations' filesep NAME_GRD_transects]; if ~exist([dirout],'dir'), mkdir(dirout); end

%% profile parameters 
h1 = 15; 
resolution_m = 1; % m 

%% READ DATA 
file_transects =['E:\West_Coast_TWL_Hazards\01_Data\Transects' filesep NAME_GRD_transects '_Transects_v3.mat'];
load(file_transects);
eval(['Transects=' NAME_GRD_transects '_Transects_v3;']);
%Define SWAN Bathymetry data file path
DirBathy='E:\West_Coast_TWL_Hazards\01_Data\SWAN_Grids\Oregon_All_grids';
Dir_Shoreline=['E:\West_Coast_TWL_Hazards\_STEP' filesep '_' NAME_GRD_transects filesep NAME_GRD_transects '_Shoreline_50m_Smooth.mat'];
load(Dir_Shoreline);
eval(['Shoreline=' NAME_GRD_transects '_Shoreline_50m_Smooth;']);
Shoreline=[Shoreline.X;Shoreline.Y];

NAME_GRD_Bathy       ={'LDUTM_100m_1','LDUTM_100m_2'};


%% FIND SWAN GRID FOR EACH TRANSECT 
SWANgrdcenter=zeros(numel(NAME_GRD_Bathy),2); 
disp('reading swan grids...') 
clear SWAN GRID SWANgrd*

for grid =1:numel(NAME_GRD_Bathy)
    
    FILENAME = [DirBathy filesep NAME_GRD_Bathy{grid} '.grd'];
    
    %wlgrid is from the Deltares Matlab toolbox
    [X,Y,ENC] = wlgrid('read', FILENAME);
       
    SWANgrdcenter(grid,1)=mean([min(X(:)) max(X(:))]); 
    SWANgrdcenter(grid,2)=mean([min(Y(:)) max(Y(:))]); 
    
    SWAN{grid}.X= X;  
    SWAN{grid}.Y= Y;     
end

%% find which grid to use 
disp('finding right grids...') 

coords = [Transects.X(:,2) Transects.Y(:,2)];
Ntr     = numel(coords(:,1));  % update Ntr to number of head of profiles that are not nan 

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


% FIND CORRECT GRID 
GRIDid = nan(Ntr,1); 
disp('finding right grid for each point...') 
for ii =1:Ntr
    VALUE =nan(numel(NAME_GRD_Bathy),1); row=VALUE; col=VALUE; 
    for grid = 1:numel(NAME_GRD_Bathy)
        if INPOINTS(ii,grid)==0; 
            continue, 
        else
            [VALUE(grid), row(grid),col(grid)] = nearestpoint_coords(SWAN{grid}.X,SWAN{grid}.Y, coords(ii,1),coords(ii,2)); 
        end % outside the grid
        
    end
    if 1-all(isnan(VALUE))
        [temp,GRIDid(ii)]=nanmax(VALUE); 
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


%Load the Bathymetries based on the grid ID's

%Control for NaNs that may have some issue being found in the tight CA
%grids however, MAKE SURE THAT ALL OF THE POINTS REASONABLY FIT INTO THIER GRID
%CATEGORIES

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



    

for grid =1:numel(GridL);

    FILENAME_grid = [DirBathy filesep NAME_GRD_Bathy{GridL(grid)} '.grd'];
    FILENAME_bathy = [DirBathy filesep NAME_GRD_Bathy{GridL(grid)} '.dep'];

    Gri = wlgrid('read', FILENAME_grid);
    Bath = wldep('read', FILENAME_bathy,Gri);
    
    NewGroup.X=Gri.X;
    NewGroup.Y=Gri.Y;
    NewGroup.Bathy=Bath(1:end-1,1:end-1);
    
    %convert to xyz to use scattered interpolant
    xx=[];
    yy=[];
    zz=[];
    for mm=1:length(NewGroup.X(1,:));
        for nn=1:length(NewGroup.X(:,1));
            xt=NewGroup.X(nn,mm);
            xx=[xx;xt];
            yt=NewGroup.Y(nn,mm);
            yy=[yy;yt];
            zt=NewGroup.Bathy(nn,mm);
            zz=[zz;zt];
            
        end
    end 
     
    NewGroup.xx=xx;
    NewGroup.yy=yy;
    NewGroup.zz=zz;
    

    eval(['Bathy_' NAME_GRD_Bathy{GridL(grid)} '=NewGroup;']);


end

%% Interpolate each transect by its index



Xp=[Transects.X];
Yp=[Transects.Y];



%% Interpolate DEM at profiles 
depth_lim = [0 h1]; % profile limits, positive is bathymetry (because pulling from DEFLT3D .dep file)
array_for_plots = [1:Ntr]; 

time1 = now; 

transects_mod = Transects; 

INCORRECT_RUNS = [];
 Bathy_name=['Bathy_' NAME_GRD_Bathy{GRIDid(1)}];
for ii = 1:Ntr%1:Ntr 
    
    OBJECTID = ii;
    disp([num2str(ii),'/',num2str(Ntr)]) 
    
    %load points separately

    Lx = Xp(ii,1)-Xp(ii,2); 
    Ly = Yp(ii,1)-Yp(ii,2); 
    
    L0 = sqrt(Lx.^2 + Ly.^2);  % true distance, from UTMs 
        
    Npts = L0./resolution_m; 
    Npts = round(Npts); 
    
    xpf   = linspace(Xp(ii,1) ,Xp(ii,2) ,Npts);
    ypf   = linspace(Yp(ii,1) ,Yp(ii,2) ,Npts);
    
    %add a bit to run an intersection and limit 
    L1=[xpf(1:end-5);ypf(1:end-5)];
    L1=[L1(1,1) , L1(1,end-5); L1(2,1) , L1(2,end-5)];
    P = InterX(Shoreline,L1);
    
    if isempty(P)
        signal=0;
    else
        signal=1;
        
        %calculate distance from these points to endpoint and grab
        %whichever is closer
        
        xcomp=(P(1,:)-xpf(end)).^2;
        ycomp=(P(2,:)-ypf(end)).^2;
        d_comp=sqrt(xcomp+ycomp);
        [~,rose]=min(d_comp);
                
        %now cut the profiles' length based on the closest point in xpf,ypf
        [~,test]=min(abs(xpf-P(1,rose)));
        test=test+1;
        xpf=xpf(test:end);
        ypf=ypf(test:end);
    end
    
   
    if strcmp(Bathy_name,['Bathy_' NAME_GRD_Bathy{GRIDid(ii)}])==0 || ii==1;
        eval(['quick_bath=Bathy_' NAME_GRD_Bathy{GRIDid(ii)} ';']);
        XX=quick_bath.xx;
        YY=quick_bath.yy;
        ZZ=quick_bath.zz;
    
        tempx=find(isnan(XX));
     
        ZZ(tempx)=[];
        XX(tempx)=[];
        YY(tempx)=[];
    
        tempz=find(ZZ==-999);
        ZZ(tempz)=[];
        XX(tempz)=[];
        YY(tempz)=[];
    
        %Scattered Interpolant method to get the nearest points
        F = scatteredInterpolant(XX,YY,ZZ,'nearest');
    end
    depthi = F(xpf,ypf);   
    dep_end=numel(depthi);
    
    
    %New addition: If the transect does not make it to land (or some other depth) for some
    %reason, extend the transect to the onshore point;
    
    if mean(depthi)==-9 && signal==0;
        %need to extend transect offshore
        D=sqrt(((xpf(1)-xpf(end))^2)+((xpf(1)-xpf(end))^2));
        
        
        U=xpf(1)-xpf(end);
        V=ypf(1)-ypf(end);
        [r,az] = pcoord(U,V);
        r=r+3000;
        [u,v] = rcoord(r, az);
        
        X1=xpf(1)+u;
        Y1=ypf(1)+v;
        
        
        Lx = X1-xpf(end); 
        Ly = Y1-ypf(end); 
    
        L0 = sqrt(Lx.^2 + Ly.^2);  % true distance, from UTMs 
        
        Npts = L0./resolution_m; 
        Npts = round(Npts); 
        
        xpf   = linspace(X1 ,xpf(end) ,Npts);
        ypf   = linspace(Y1 ,ypf(end) ,Npts);
        depthi = F(xpf,ypf);
    end
        
        
    
    if depthi(end)>10; %in meters

        Lx = Xp(ii,1)-Xp(ii,3); 
        Ly = Yp(ii,1)-Yp(ii,3); 
    
        L0 = sqrt(Lx.^2 + Ly.^2);  % true distance, from UTMs 
        
        Npts = L0./resolution_m; 
        Npts = round(Npts); 
    
        xpf   = linspace(Xp(ii,1) ,Xp(ii,3) ,Npts);
        ypf   = linspace(Yp(ii,1) ,Yp(ii,3) ,Npts);
        depthi = F(xpf,ypf);
        
        if depthi(end)<10; %in meters
            test_depth=depthi(dep_end:end);
            a=find(test_depth<0);
            if ~isempty(a);
            a=a(1);
            depthi((dep_end+a):end)=[];
            end
        end
        
    end
        
        
   
%%% Now pick out the index of the seaward most
%%% point at 15m depth. The issue here will be if it crosses something
%%% that's shallow or if its in a bay. 
    if signal==1;
        I = GetContourPoint_OR(depthi,h1,50,resolution_m);
    else
       I = GetContourPoint_OR_no_rocks(depthi,h1,50,resolution_m);
    end
    
    if depthi(I)<0;
        I=I+1;
        if depthi(I)<0;
            I=I-2;
        end
    end
        
    
    
    try 
        I2=I-1;
        
        Ii=abs(h1-depthi(I));
        I2i=abs(h1-depthi(I2));
        
        if I2i<Ii
            I=I2;
        end
    end
        
       
    
    

    
    if I==0;
    Xpu=NaN;
    Ypu=NaN;
    BathyContour.X(ii)=NaN;
    BathyContour.Y(ii)=NaN;   
    BathyContour.Z(ii)=NaN;    
        
    else    
    Xpu=xpf(I);
    Ypu=ypf(I);
    BathyContour.X(ii)=Xpu;
    BathyContour.Y(ii)=Ypu;   
    BathyContour.Z(ii)=depthi(I);
    end
    Bathy_name=['Bathy_' NAME_GRD_Bathy{GRIDid(ii)}];
   
    end

    %% SAVE  
    figure;scatter(BathyContour.X, BathyContour.Y);
    name=[NAME_GRD_transects '_Contour_pts'];
    eval([name '=BathyContour;']);
    fout = [dirout, filesep, name '.mat']; 
    save(fout,name);
    
    

    



