%%Third Scipt to calculate profile elevations
warning off all 
%% CONFIG 
NAME_GRD_transects  = 'Douglas'; 
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

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
% read transects  

%%
num_dem=1;
%%

load('Transects2.mat');
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


% load R-Topo and Fbount
for rr=1:num_dem;
    load(['Rtopo_' num2str(rr)]);
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




% 2) FIND RIGHT GRID 

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

close all
GridL=unique(GRIDid);



depth_lim = [h1 h2]; % profile limits, bathymetry <0 

array_for_plots = [0]; 

time1 = now; 

transects_mod = Transects2;

INCORRECT_RUNS = [];

for ii = 1:Ntr;
    
    if isnan(transects_mod.xpf(ii,1));
        continue
    end
    OBJECTID = transects.OBJECTID(ii); 
    disp([num2str(ii),'/',num2str(Ntr)]) 
    Lx = transects_mod.xpf(ii,1)-transects_mod.xpf(ii,2); 
    Ly = transects_mod.ypf(ii,1)-transects_mod.ypf(ii,2); 
    
    L0 = sqrt(Lx.^2 + Ly.^2);  % true distance, from utms 
        
    Npts = L0./resolution_m; 
    Npts = round(Npts); 
    
    %get the first intersection from X0,Y0 to the boundary
    
    
    xpf   = linspace(transects_mod.xpf(ii,1) ,transects_mod.xpf(ii,2) ,Npts);
    ypf   = linspace(transects_mod.ypf(ii,1) ,transects_mod.ypf(ii,2) ,Npts);

    uu=GRIDid(ii);
    Int_Name=['F_' num2str(uu)];
    if ~exist(Int_Name,'var')
        load(Int_Name)
        
        for jj=1:num_dem;
            if jj~=uu;
            Bad_Name=['F_' num2str(jj)];
            eval(['clear ' Bad_Name]);
            end
        end
    end
    
    eval(['depthi = F_' num2str(uu) '(xpf,ypf);']);
   
    
    if sum(~isnan(depthi))==0 || isempty(depthi);
        continue
    end

    indland=depthi>=MHW(ii,3); 
    indN=find(~isnan(depthi));
         
    %separate the results by distance
    
    diffinland=diff(indland); 
    %define blocks of land 
    ind1 = find(diffinland==1); 
    ind_1= find(diffinland==-1); 
    
    
    if numel(ind1)>numel(ind_1);
        a=numel(ind1)-numel(ind_1);
        ind1(end-a+1:end)=[];
    elseif numel(ind1)<numel(ind_1);
        a=numel(ind_1)-numel(ind1);
        ind_1(end-a+1:end)=[];
    end
    
    if numel(ind1)>1% more than 1 intersection with land
        % 1 - check length of island
      
        distances = (ind_1-ind1)*resolution_m;
        % 2 - check if close to initial contour 
%  
        
        for yy=1:numel(ind1)
            if distances(yy)<prof_cut_len
                if yy==numel(ind1)
                else
                depthi(ind1(yy)-1:ind_1(yy)+1)=1;
                end
            else 
                break %this bit should keep everything past the inital large step from being removed
            end
        end
        

    end
        
        
        
    coastal = depthi>=depth_lim(1) & depthi<=depth_lim(2);
    
    if max(depthi(coastal)) < depth_lim(2); % CORRECT FOR LOWER ELEVATIONS
        [IND, D] = nearestpoint(depth_lim(2),depthi);
        coastal(IND)=1; 
    end
    
    if any(ii == array_for_plots)  
        figure('visible',VISIBILITY) 
        plot(xpf,depthi) 
        hold on, plot(xpf(coastal), depthi(coastal),'.r') 
        axis tight 
        
    end      
        
    %% Edit transect for non hydroflattened data     
    depth=depthi;
  
   
    %Remove any potential islands in front of the profile by finding any
    %regions that are above MSL but smaller than 50m 
   
    temp_0=find(depth<MSL(ii,3));
    if ~isempty(temp_0);
    temp=diff(temp_0);
    test=find(temp>50)-1;
    
    qwe=find(test==0);
    if ~isempty(qwe)
        test(qwe)=1;
    end
    
    if~isempty(test)
    temp(test:end)=1;
    end
    temp(temp>50)=1;
    temp=find(temp>3);
    
    
    temp=temp_0(temp);
    if ~isempty(temp)
    depth(1:temp(end))=NaN;
    else
        if ~isempty(test)
    depth(1:temp_0(test(1)))=NaN; 
        end
   
    end
    end
    
   %Want everything seaward of some base height (~3 m)
   IndMax=find(depth>3);
   if ~isempty(IndMax)
   IndMax=IndMax(1); %Grab the first point that is greater than the threshold
   else 
       try
        IndMax=find(depth>2);
        IndMax=IndMax(1);
       catch
           try
         IndMax=find(depth>0.3);  
         IndMax=IndMax(1);
           catch
               disp('All Water');
               continue
           end
       end
   end
    
   cutoff =MHW(ii,3); %Anything above this mark will not be considered 
   cutidx=find(depth>cutoff);
   testdepth=depth;
    testdepth(IndMax:end)=NaN;
   testdepth(cutidx)=NaN;
    
    
    d=diff(testdepth(1:2:end));
    didx=max(find(d<0))*2;
    
   
    xpf(1:didx)=NaN;
    ypf(1:didx)=NaN;
    depth(1:didx)=NaN;
    
    
    
    
    %%    
    if isempty(depth), continue, end    
        
    L = 0:resolution_m:(numel(depth)-1).*resolution_m; 
   
    %% ISLANDS FILTER 
    hend = depth(end); 
    if hend<0 && any(depth(:)>0)
        ind0=find((depth)>=0); 
        depth(ind0(end):end)=[]; 
        %bedtype(ind0(end):end)=[]; 
        xpf(ind0(end):end)=[]; 
        ypf(ind0(end):end)=[]; 
        L(ind0(end):end)=[]; 
    end
        
    check_plot = 0;    
        
     %% NO LAND FLAG
    if all(depth<0), disp('NO LAND'), continue, end 
    %% figures     
    if any(ii == array_for_plots)  
%         close all 
        figure, plot(xpf,depth) 
        set(gca,'ylim',[-50 50])
        %grid on 
        hold on, plot(xpf, depth,'.k') 
        scatter(xpf,depth,20,'filled')
        
        legend('DEM','selected profile','active profile')
        legend('location','southeast') 
        
        plot(xpf,xpf.*0,'-K','linewidth',2), text(xpf(1) ,5,'MSL') 
        title(['Profile ',dec2base(ii,10,4),' - Grid: ',NAME_GRD_transects])

    end
    if rem(ii,5) == 0, close all, end 
    
    %% save profile data
    wave_profile=[]; 
    
    
    %remove NaN values
    w=hanning(2*resolution_m+1);
    depth=filter(w/sum(w),1,depth);
    lag=((2*resolution_m+1-1)/2);
    depth=depth(1+lag:end);    

    indN=find(~isnan(depth));
    
    wave_profile.depth    = depth(indN);
    
    
    xpf=xpf(1:end-lag);
    wave_profile.xpf      = xpf(indN); 
    ypf=ypf(1:end-lag);
    wave_profile.ypf      = ypf(indN); 
    wave_profile.true_res = resolution_m; 
    L=L(1:end-lag);
    wave_profile.L        = L(indN); 
    wave_profile.OBJECTID = OBJECTID;
    
    
    wave_profile.Lcoast = wave_profile.L(1);
    wave_profile.xcoast = wave_profile.xpf(1); 
    wave_profile.ycoast = wave_profile.ypf(1); 

%add MHW and MSL elevations to data structure
     mhw=cutoff;

    
    temp=find(wave_profile.depth>mhw);
    if ~isempty(temp);
    temp=[temp(1)-1 temp(1)];
    if temp(1)<1
        I=1;
    else
    val=wave_profile.depth(temp);
    val=abs(val-mhw);
    [m I]=min(val);
    end
    
    wave_profile.mhw=mhw;
    wave_profile.mhw_Lpos=temp(I);
    else
     wave_profile.mhw=mhw;
    wave_profile.mhw_Lpos=wave_profile.L(end);   
    end
    
 msl=MSL(ii,3);    
    temp=find(wave_profile.depth>msl);
    if ~isempty(temp);
    temp=[temp(1)-1 temp(1)];
    if temp(1)<1
        I=1;
    else
    val=wave_profile.depth(temp);
    val=abs(val-msl);
    [m I]=min(val);   
    end
    wave_profile.msl=msl;
    wave_profile.msl_Lpos=temp(I);
    else
        wave_profile.msl=msl;
        wave_profile.msl_Lpos=wave_profile.L(end);
    end


    %% save 
    fout = [dirout, filesep, 'profile_', dec2base(OBJECTID,10,4) ,'.mat']; 
    m=matfile(sprintf('%s', fout),'writable',true); 
    m.wave_profile = wave_profile; 


  
    
end
    
    winopen(dirout) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
