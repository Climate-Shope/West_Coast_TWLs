%% Script to find erroneous runup methods saved to profile files and 
%  correct them 
%% Define parameters
clear; close all; run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = {'San_Diego','Orange','Los_Angeles','Ventura','Santa_Barbara',...
    'Santa_Barbara_West','San_Luis_Obispo','Monterey','Santa_Cruz','San_Mateo','San_Francisco',...
    'Marin','Sonoma','Mendocino','Humboldt','Del_Norte','Curry','Coos','Douglas','Lane',...
    'Lincoln','Tillamook','Clatsop','Pacific','Grays_Harbor','Jefferson','Clallam'}; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

UTM={'11 N','11 N','11 N','11 N','11 S','10 S','10 N','10 N','10 N',...
     '10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N',...
     '10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N'};
 
Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);
 
%% Loop through each location
 for WestCoast=1:27
     load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\time_series_reconstructed',filesep,'ii_grid_array.mat']);
    Ntr=ii_grid_array;
    
    DirOut=['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'Morphology_8_29_2019'];
    if ~exist(DirOut,'dir'), mkdir(DirOut); end
    %Load in correct ESI value 
    load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\Profile_ESI',filesep,NAME_GRD_transects{WestCoast}, '_ESI_Array.mat']);
     cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');

%% Loop though each transect  
     for ii=ii_grid_array;
    OBJECTID=dec2base(ii,10,4);
    
    if exist(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'Morphology_8_20_2019', filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat'],'file');
    load(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'Morphology_8_20_2019', filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat']);
    Rey=SlopeInfo.Rey;Flag for seawall
    %If profile is poorly interpolated, skip it 
    if SlopeInfo.maxi==0;
        continue
    end
    %% Load in base Profile
    load(['F:\West_Coast_TWL_Hazards\03_Results\',NAME_GRD_transects{WestCoast},filesep, 'interp_profiles_wave_model_edited', filesep, 'profile_', dec2base(ii,10,4), '.mat']);
    
    %Santa Barabara transect list is revered
    if strcmp(NAME_GRD_transects{WestCoast},'Santa_Barbara')
       wave_profile.depth=wave_profile.depth'; 
    end
    
    %clip to 300m onshore
    if wave_profile.msl_Lpos>300
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
    a=find(wave_profile.depth>30);
    if ~isempty(a);
        a=a(1);
        
        if a>10
        wave_profile.depth=wave_profile.depth(1:a);
        wave_profile.xpf=wave_profile.xpf(1:a);
        wave_profile.ypf=wave_profile.ypf(1:a);
        wave_profile.L=wave_profile.L(1:a);
        end
    end
    
    ESI=SlopeInfo.ESI;

%% Simplify profile 
     result = DouglasPeucker([wave_profile.L; wave_profile.depth],1.0);
     if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4')
    if strcmp(ESI,'6B') && ~strcmp(ESI_Array(str2num(OBJECTID),5),'3A: Fine- to medium-grained sand beaches');
    [pks,locs] = findpeaks(result(2,:),0.2,1,3.5);
    
    elseif strcmp(ESI,'6B') && strcmp(ESI_Array(str2num(OBJECTID),5),'3A: Fine- to medium-grained sand beaches');
    [pks,locs] = findpeaks(result(2,:),0.3,1,3.5);
    else
    [pks,locs] = findpeaks(result(2,:),1,1,3.5); %Must have a 1 m prominence 
    end
    depth=interp1(result(1,:),result(2,:),wave_profile.L);  
    if ~isempty(locs)
    locs=locs(1);
    locs=find(wave_profile.L==result(1,locs));
    end
    
    if ~isempty(locs)
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
    end
        
        
        
     end
   
   %%%% 4/26/2019 Added  
if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4')
         avg=wave_profile.mhw;
       L1=[[1:length(wave_profile.depth)];wave_profile.depth];
       L2=[[1 length(wave_profile.depth)];[avg avg]];
        P = InterX(L1,L2);
        if ~isempty(P);
            T=round(P(1,:)); 
            T(T>30)=[];
            if ~isempty(T);
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



 %% Pick the correct average Beach Slope for the profile, doesnt matter what the previous was determiend to be
   %first determine the indicies where the correct county name 
    uwu=find(strcmp(AllStockdon.Counties,NAME_GRD_transects{WestCoast}));
    %second determine the indicies where the correct profile name
    uvu=find(strcmp(AllStockdon.Names,['profile_' dec2base(ii,10,4)]));
    %now find the overlap,  should really be one value
    uzu=find(ismember(uwu,uvu));
    uzu=uwu(uzu);
    fslope=AllStockdon.Slope_MeanMajors_Limited(uzu);
    fslope=repmat(fslope,1,3); %repmat for compatibility with older functions 


 %% Update the runup method 
 
    [Runup_Method]=Runup_Method_Update(SlopeInfo,wave_profile,fslope);
    Runup_Method.StockSlope=fslope(1);
    SlopeInfo.Runup_Method=Runup_Method;
    
    if ~isempty(SlopeInfo.Rey);
        if SlopeInfo.toe_loc<1
        slope=(SlopeInfo.depth(SlopeInfo.maxi)-SlopeInfo.depth(1))/(SlopeInfo.maxi-1);
        slope=atan(slope)*(180/pi);
        else
        
        slope=(SlopeInfo.depth(SlopeInfo.maxi)-SlopeInfo.depth(SlopeInfo.toe_loc))/(SlopeInfo.maxi-SlopeInfo.toe_loc);
        slope=atan(slope)*(180/pi);
        end
        if slope<10 && strcmp(SlopeInfo.Runup_Method.Method,'Stockdon') && (SlopeInfo.depth(SlopeInfo.maxi)-SlopeInfo.depth(SlopeInfo.toe_loc))<1.5;
            SlopeInfo.Runup_Method.ExceedMethod='Stockdon';
        end
        
    end
    
    name=[DirOut filesep 'SlopeInfo_Transect_' num2str(ii)]; 
    save(name,'SlopeInfo');
    
    clear SlopeInfo wave_profile
    end
    end
 end