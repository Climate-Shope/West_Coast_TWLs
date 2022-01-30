%Code to automatically extract relevant elevation information from elevation profiles
%Code listed below is region specific and may require editing
%for another region

%% Load directories 
clear; close all; run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Douglas'; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI', filesep, NAME_GRD_transects, '_ShoreType.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_SeaIntersect.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_ESI_Array.mat']);

load([directories.timeseries filesep 'ii_grid_array']);


utmzone='10 N'; %Needed for lat lon transformations 10N for Oregon

%Load in the Precalculated Stockdon-type Slopes calculated separately
Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);
% Create the out directory 

DirOut=['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Morphology'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end

cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');%Where these functions are housed
for ii=ii_grid_array;

    OBJECTID=dec2base(ii,10,4);
    
    %Do not need wave conditions  
       
    % load the profile 
    try
    load([directories.profiles, '_edited', filesep, 'profile_', dec2base(ii,10,4), '.mat']);
    catch
     continue   
    end
    
    if length(wave_profile.depth(:,1))>1
        wave_profile.depth=wave_profile.depth';
    end
    
     %Skip the loop if the profile is empty
    if isempty(wave_profile.depth)
        SlopeInfo.maxi=0;
        SlopeInfo.maxi_ele=0;
        SlopeInfo.MHW=wave_profile.mhw;
        SlopeInfo.MHW_loc=0;
        SlopeInfo.MSL=wave_profile.msl;
        SlopeInfo.CliffTop_Ele=0;
        SlopeInfo.CliffTop_Loc=0;
        SlopeInfo.CliffJun_Loc=0;
        SlopeInfo.CliffJun_Ele=0;
        SlopeInfo.toe_loc=0; 
        SlopeInfo.toe_ele=0;
        SlopeInfo.toe_onshore=0; 
        SlopeInfo.maxi_onshore=0; 
        SlopeInfo.berm_loc=0;
        SlopeInfo.berm_width=0;
        SlopeInfo.berm_toe=0;
        SlopeInfo.ESI=ESI;
        SlopeInfo.extra_toe=0;
        SlopeInfo.Rey=Rey;
        SlopeInfo.depth=wave_profile.depth;
        SlopeInfo.xpf=wave_profile.xpf;
        SlopeInfo.ypf=wave_profile.ypf;
        SlopeInfo.L=wave_profile.L;
        SlopeInfo.overtop_point=0;
        name=[DirOut filesep 'SlopeInfo_Transect_' num2str(ii)]; 
        save(name,'SlopeInfo');
    
        clear SlopeInfo wave_profile
        continue
    end

    
  
    
    if wave_profile.msl_Lpos>300 %Cut to just 300m onshore
        wave_profile.depth=wave_profile.depth(wave_profile.msl_Lpos-1:end);
        wave_profile.xpf=wave_profile.xpf(wave_profile.msl_Lpos-1:end);
        wave_profile.ypf=wave_profile.ypf(wave_profile.msl_Lpos-1:end);
        wave_profile.L=wave_profile.L(wave_profile.msl_Lpos-1:end);
        wave_profile.Lcoast= wave_profile.Lcoast-wave_profile.msl_Lpos-1;
        wave_profile.mhw_Lpos=wave_profile.mhw_Lpos-wave_profile.msl_Lpos+2;
        wave_profile.msl_Lpos=wave_profile.msl_Lpos-wave_profile.msl_Lpos+2;
    end
    
     if length(wave_profile.depth)>300 && max(wave_profile.depth(1:300))>10
        wave_profile.depth=wave_profile.depth(1:300);
        wave_profile.xpf=wave_profile.xpf(1:300);
        wave_profile.ypf=wave_profile.ypf(1:300);
        wave_profile.L=wave_profile.L(1:300);

    end
    
    %Cut the length based on the onshore height
    a=find(wave_profile.depth>30); %Anything greater than 30m is in no danger of overtopping
    if ~isempty(a);
        a=a(1);
        
        if a>10
        wave_profile.depth=wave_profile.depth(1:a);
        wave_profile.xpf=wave_profile.xpf(1:a);
        wave_profile.ypf=wave_profile.ypf(1:a);
        wave_profile.L=wave_profile.L(1:a);
        end
    end
    
    
  %Skip the loop if the profile is empty
    if isempty(wave_profile.depth)
        SlopeInfo.maxi=0;
        SlopeInfo.maxi_ele=0;
        SlopeInfo.MHW=wave_profile.mhw;
        SlopeInfo.MHW_loc=0;
        SlopeInfo.MSL=wave_profile.msl;
        SlopeInfo.CliffTop_Ele=0;
        SlopeInfo.CliffTop_Loc=0;
        SlopeInfo.CliffJun_Loc=0;
        SlopeInfo.CliffJun_Ele=0;
        SlopeInfo.toe_loc=0; 
        SlopeInfo.toe_ele=0;
        SlopeInfo.toe_onshore=0; 
        SlopeInfo.maxi_onshore=0; 
        SlopeInfo.berm_loc=0;
        SlopeInfo.berm_width=0;
        SlopeInfo.berm_toe=0;
        SlopeInfo.ESI=ESI;
        SlopeInfo.extra_toe=0;
        SlopeInfo.Rey=Rey;
        SlopeInfo.depth=wave_profile.depth;
        SlopeInfo.xpf=wave_profile.xpf;
        SlopeInfo.ypf=wave_profile.ypf;
        SlopeInfo.L=wave_profile.L;
        SlopeInfo.overtop_point=0;
        name=[DirOut filesep 'SlopeInfo_Transect_' num2str(ii)]; 
        save(name,'SlopeInfo');
    
        clear SlopeInfo wave_profile
        continue
    end
    
    
    
    
    if ii<length(ShoreType.ESI)
       ESI=ShoreType.ESI{ii};
    else
       ESI=ShoreType.ESI{end}; 
    end
    
    
    
    
    
    
    %% Simplify Profile   
    %These simplification values are calibrated for each region and ESI, but the
    %user can decide the appropriate threshold
    
    
    
    
    %Simplifying profile using DouglasPeucker should not be as strong for a
    %low-lying profile
    if max( wave_profile.depth) < 5.5
        result = DouglasPeucker([wave_profile.L; wave_profile.depth],0.1);
    else
        result = DouglasPeucker([wave_profile.L; wave_profile.depth],1.0);
    end
    
    %% Correct profile length to account for engineered and rocky shorelines
    
    
    if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4') %not sand or gravel beach
    
    %Cut the profile to the correct length by defining an onshore location
    if strcmp(ESI,'6B') && ~strcmp(ESI_Array(str2num(OBJECTID),5),'3A: Fine- to medium-grained sand beaches'); %Just riprap
        [pks,locs] = findpeaks(result(2,:),0.2,1,3.5);
    
    elseif strcmp(ESI,'6B') && strcmp(ESI_Array(str2num(OBJECTID),5),'3A: Fine- to medium-grained sand beaches'); %riprap with a sandy beach
        [pks,locs] = findpeaks(result(2,:),0.3,1,3.5);
    
    %Again, lower-lying profiles have slightly different parameters 
    elseif max( wave_profile.depth) < 5.5
        [pks,locs] = findpeaks(result(2,:),[],1,3.5); %Must have a 1 m prominence
    else
        [pks,locs] = findpeaks(result(2,:),1,1,3.5); %Must have a 1 m prominence
    end
    
    %Iterpolate the simplified profile to the same number of points as the
    %original elevation profile
    depth=interp1(result(1,:),result(2,:),wave_profile.L);  
    
    if ~isempty(locs) %If a shoreward maxima is found 
        
        if numel(locs)>1
            first=find(wave_profile.L==result(1,locs(1)));
            second=find(wave_profile.L==result(1,locs(2)));
            if second-first>200 %If its really far apart, more likely to be a harbor
                careful=1;
            end
        end
    locs=locs(1);
    locs=find(wave_profile.L==result(1,locs));
    end
    
    %now if there are maxima and a harbor situaiton is not detected
    if ~isempty(locs) && ~exist('careful','var')
        %account for low lying conditions
    if locs(1)<10 && wave_profile.depth(locs(1))<5 
        ee=locs(1);
        a=diff(depth);a=diff(a);
        [pks,locs] = findpeaks(a,[],1,0.1);
        locs(locs<ee)=[];
        if ~isempty(locs)
        cut=locs(1)+1;
        wave_profile.depth=wave_profile.depth(cut:end);
        wave_profile.xpf=wave_profile.xpf(cut:end);
        wave_profile.ypf=wave_profile.ypf(cut:end);
        wave_profile.L=wave_profile.L(1:length(wave_profile.depth));
        end
    end
    else
       
    end
        
        
        
    end
   
   %%%% Additional consideration if the profile does not intersect with
   %%%% calculated MHW due to missing data (also happens along plunging
   %%%% cliffs)
if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4')
       avg=wave_profile.mhw;
       L1=[[1:length(wave_profile.depth)];wave_profile.depth];
       L2=[[1 length(wave_profile.depth)];[avg avg]];
       P = InterX(L1,L2);
        if ~isempty(P);
            T=round(P(1,:)); 
            T(T>30)=[];
            if ~isempty(T) && ~exist('careful','var');
            if length(T)>1;
                T=T(end);
            wave_profile.depth(1:T-1)=[];
            wave_profile.xpf(1:T-1)=[];
            wave_profile.ypf(1:T-1)=[];
            wave_profile.L(1:T-1)=[];
            end
            end
        end
end
     %%%%
     
 if exist('careful','var');   
     clear careful
 end
 
 %% Fill in the ESI shoreline type and seawall presence for the combined structure file
    ShoreID=find(ShoreType.OBJECTID==ii);
    if isempty(ShoreID);
        ESIID={NaN};
    elseif isnan(ShoreType.ESI{ShoreID,1})
        ESIID={NaN};
    else
        ESIID=cellstr(ShoreType.ESI{ShoreID,1});
    end
    
    %check for a detached seawall
    if ~isempty(SeaIntersect);
    temp=find(SeaIntersect(:,1)==ii);
    %Rey variable carries across subsequent files to represent a detached
    %seawall
    if ~isempty(temp);
        Rey=SeaIntersect(temp,:);
    else
        Rey=[];
    end
    else
        Rey=[];
    end
    if ii<length(ShoreType.ESI)
     ESI=ShoreType.ESI{ii};
    else
     ESI=ShoreType.ESI{end}; 
    end
 

    %% Pick the correct average Stockdon Slope for the profile (calculated externally), doesnt matter what the previous was determiend to be
   %first determine the indicies where the correct county name 
    uwu=find(strcmp(AllStockdon.Counties,NAME_GRD_transects));
    %second determine the indicies of the correct profile name
    uvu=find(strcmp(AllStockdon.Names,['profile_' dec2base(ii,10,4)]));
    %now find the overlap, should normally on be position one 
    uzu=find(ismember(uwu,uvu));
    uzu=uwu(uzu);
    fslope=AllStockdon.Slope_MeanMajors_Limited(uzu);
    fslope=repmat(fslope,1,3);
    
   
   %% Skip the loop if the profile maximum does not exceed MHW
   tt=sort(wave_profile.depth);pct=round(0.5*length(wave_profile.depth));
    if mean(tt(pct:end))<wave_profile.mhw
        SlopeInfo.maxi=0;
        SlopeInfo.maxi_ele=0;
        SlopeInfo.MHW=wave_profile.mhw;
        SlopeInfo.MHW_loc=0;
        SlopeInfo.MSL=wave_profile.msl;
        SlopeInfo.CliffTop_Ele=0;
        SlopeInfo.CliffTop_Loc=0;
        SlopeInfo.CliffJun_Loc=0;
        SlopeInfo.CliffJun_Ele=0;
        SlopeInfo.toe_loc=0; 
        SlopeInfo.toe_ele=0;
        SlopeInfo.toe_onshore=0; 
        SlopeInfo.maxi_onshore=0; 
        SlopeInfo.berm_loc=0;
        SlopeInfo.berm_width=0;
        SlopeInfo.berm_toe=0;
        SlopeInfo.ESI=ESI;
        SlopeInfo.extra_toe=0;
        SlopeInfo.Rey=Rey;
        SlopeInfo.depth=wave_profile.depth;
        SlopeInfo.xpf=wave_profile.xpf;
        SlopeInfo.ypf=wave_profile.ypf;
        SlopeInfo.L=wave_profile.L;
        SlopeInfo.overtop_point=0;
        name=[DirOut filesep 'SlopeInfo_Transect_' num2str(ii)]; 
        save(name,'SlopeInfo');
    
        clear SlopeInfo wave_profile
        continue
    end
   
    
    %% provide a fix where the seawall is very far removed from the profile
    if ~isempty(Rey)
    if numel(num2str(round(Rey(:,3))))<3
      xx= wave_profile.xpf;
      yy= wave_profile.ypf;
      utmzone=repmat(wave_profile.UTMZONE,length(xx),1); 
       [Lat,Lon] = utm2deg(xx,yy,utmzone);
       xcomp=(Rey(:,2)-Lon).^2;
       ycomp=(Rey(:,3)-Lat).^2;
    else
       xcomp=(Rey(:,2)-wave_profile.xpf).^2;
       ycomp=(Rey(:,3)-wave_profile.ypf).^2;   
    end
   
   

    distance=sqrt(xcomp+ycomp);
    val=min(min(distance));

    %If far away, could be wrong location or irrelevant to calculations
    if val>100
        Rey=[];
    end
    end
    
    
    
    
    
    
   
   %% extract profile parameters, calibrated by region
    [SlopeInfo]=MorphLove_Oregon(wave_profile,ESI,ESI_Array,Rey,OBJECTID);

  %Provide a reanalysis if the analyzed morphology produced a toe location
  %further onshore than the slopebreak of a cliff and recut profile length
  %to be based on the locaiton of the MHW intersection
    if SlopeInfo.toe_ele > SlopeInfo.CliffJun_Ele; %If this is the case, find the next MHW positions
       L1=[[1:length(wave_profile.depth)];wave_profile.depth];
       L2=[[1 length(wave_profile.depth)];[wave_profile.mhw wave_profile.mhw]];
        P = InterX(L1,L2);
        if ~isempty(P);
            P=round(P);
            try
                P=P(1,2);
            catch
                P=P(1,1);
            end
            wave_profile.depth=wave_profile.depth(P:end);
            wave_profile.xpf=wave_profile.xpf(P:end);
            wave_profile.ypf=wave_profile.ypf(P:end);
            wave_profile.L=wave_profile.L(P:end);
        end
        %rerun morphology analysis code
        [SlopeInfo]=MorphLove_Oregon(wave_profile,ESI,ESI_Array,Rey,OBJECTID);
    end
        
        
    Rey=SlopeInfo.Rey; %flag and position of detached seawall 
 
%% calculate run-up methods to be used at the profile     
    
    [Runup_Method]=Runup_Method(SlopeInfo,wave_profile,Rey,fslope);
    %reset the stockdon runup slope to be the average alongshore slope for
    %the locaiton of interest
    Runup_Method.StockSlope=fslope(1);
    SlopeInfo.Runup_Method=Runup_Method;
    
    %Add in fixes
    
    %For small, close to the ground seawalls, hydrodynamics are closer to stockdon so keep with that 
    if ~isempty(SlopeInfo.Rey);
        if SlopeInfo.toe_loc<1
            rad=1;
        else
            rad=SlopeInfo.toe_loc;
        end
        slope=(SlopeInfo.depth(SlopeInfo.maxi)-SlopeInfo.depth(rad))/(SlopeInfo.maxi-rad);
        slope=atan(slope)*(180/pi);
        if slope<10 && strcmp(SlopeInfo.Runup_Method.Method,'Stockdon') && (SlopeInfo.depth(SlopeInfo.maxi)-SlopeInfo.depth(rad))<1.5;
            SlopeInfo.Runup_Method.ExceedMethod='Stockdon';
        end
        
    end
    
    %Now for stockdon runup calculations, the overtopping point
    %should be the maximum, not the barrier junction or the slope break for
    %a cliff or dune
    
    if strcmp(SlopeInfo.Runup_Method.Method,'Stockdon');
        el=(SlopeInfo.depth(SlopeInfo.maxi)-SlopeInfo.depth(SlopeInfo.CliffJun_Loc));
        dd=SlopeInfo.maxi-SlopeInfo.CliffJun_Loc;
         
        %If the maximum and slope break locations are not far apart and the
        %slope break is lower than 7m (determined by testing)
         if strcmp(SlopeInfo.Runup_Method.ExceedMethod,'Stockdon')&& dd<70 && SlopeInfo.depth(SlopeInfo.CliffJun_Loc)<7;
             SlopeInfo.overtop_point=SlopeInfo.maxi;
         %If there is a sandy beach shoreline, the slope break is less than 7m, the maximum elevaiton is < 18 m, and the maximum and slope break are not too far apart
         % or if a gravel beach, the slope break is < 7 m, the maximum is greater than 7m , they are reatively close, and the elevaiton difference between the two is not more than 10m
         %each threhsold is detemrined by testing reigonally and is highly
         %depeneent on the extracted profile morphologies. The user must
         %reevaluate these for thier location of interest       
         elseif strcmp(ESI,'3A') && SlopeInfo.CliffJun_Ele<7 && SlopeInfo.maxi_ele>7 && SlopeInfo.maxi_ele<18  && SlopeInfo.maxi-SlopeInfo.CliffJun_Loc<80|| strcmp(ESI,'4') && SlopeInfo.CliffJun_Ele<7 && SlopeInfo.maxi_ele>7 && SlopeInfo.maxi-SlopeInfo.CliffJun_Loc<100 && el<10
             SlopeInfo.overtop_point=SlopeInfo.maxi;
         end
        
    end
    
    
    name=[DirOut filesep 'SlopeInfo_Transect_' num2str(ii)]; 
    save(name,'SlopeInfo');
    
    clear SlopeInfo wave_profile
    
    
end
    
    