%Load the wave condtions for a specific profile to extract the maximum
%energy conditions
clear; close all; run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Santa_Cruz'; 
run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

%% Determine how many files to loop through
load([directories.timeseries filesep 'ii_grid_array']);
Ntr=ii_grid_array;

load(['E:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI', filesep, NAME_GRD_transects, '_ShoreType.mat']);
load(['E:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_SeaIntersect.mat']);
if ~exist([directories.slope],'dir'), mkdir(directories.slope); end

%OKay first read the files in the directory 
listing = dir(['E:\West_Coast_TWL_Hazards\03_Results\Tides_Out\',NAME_GRD_transects]);
listing=size(listing);listing=listing(1)-2; 
tide_num=linspace(1,((listing-1)*10)+1,listing);


%Load in the NTR Lookup table
load('E:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\NTR_Lookup');
%Grab the Output Coordinates
NTR_Lat=cell2mat(NTR_Lookup(:,1));
NTR_Lon=cell2mat(NTR_Lookup(:,2))*-1;

utmzone='10 N';


%load in the filler slopes

load('E:\West_Coast_TWL_Hazards\01_Data\CA_Slope_Data\Fill_Slopes.mat');


% Create the out directory 
DirOut=['E:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\StockSlopes'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end


cd('E:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');

SlopeArray=[];
for ii=ii_grid_array;
    OBJECTID=dec2base(ii,10,4);
    
    %Load wave timeseries and calculate the most energetic condtion
    Wavecon=load([directories.timeseries filesep 'WD4R_prf' dec2base(ii,10,4) '.mat']);
       
    [EI,idx]=max(Wavecon.Hs);
    test_con.Hs=Wavecon.Hs(idx);test_con.T=Wavecon.T(idx);
    
    wdepth=Wavecon.z0(1);
    if wdepth>10;
    [Ho,HoTwo]=BackCalcHo(test_con.T,wdepth,test_con.Hs);
    test_con.Hs=Ho;
    end
    
    
    
    
    %now load the profile 
    load([directories.profiles, '_edited', filesep, 'profile_', dec2base(ii,10,4), '.mat']);
    
    cut=[];

    %edit the profile in case its a harbor or there is an issue with the
    %length in front of the MSL point
    
    if wave_profile.msl_Lpos>300
        wave_profile.depth=wave_profile.depth(wave_profile.msl_Lpos-1:end);
        wave_profile.xpf=wave_profile.xpf(wave_profile.msl_Lpos-1:end);
        wave_profile.ypf=wave_profile.ypf(wave_profile.msl_Lpos-1:end);
        wave_profile.L=wave_profile.L(wave_profile.msl_Lpos-1:end);
        wave_profile.Lcoast= wave_profile.Lcoast-wave_profile.msl_Lpos-1;
        wave_profile.mhw_Lpos=wave_profile.mhw_Lpos-wave_profile.msl_Lpos+2;
        wave_profile.msl_Lpos=wave_profile.msl_Lpos-wave_profile.msl_Lpos+2;
    end
    %Grab the shoreline type
    
    ShoreID=find(ShoreType.OBJECTID==ii);
    if isempty(ShoreID);
        ESIID={NaN};
    else
        ESIID=cellstr(ShoreType.ESI{ShoreID,1});
    end
    
    %check for a detached seawall
    
    temp=find(SeaIntersect(:,1)==ii);
    if ~isempty(temp);
        Rey=SeaIntersect(temp,:);
    else
        Rey=[];
    end

    %Calc Tides
    %Find closest
    [~,temp]=min(abs(tide_num-ii));
    temp=tide_num(temp(1));
    load(['E:\West_Coast_TWL_Hazards\03_Results\Tides_Out\',NAME_GRD_transects,filesep, 'Tide_Series_',num2str(temp)]);
    tide=Tide_Series.WL;
    clear Tide_Series
    
    %Find and Load the appropriate NTR Dataset
    
    %Load the first point in the profile 
    x=wave_profile.xpf(1);
    y=wave_profile.ypf(1);
    
    [Lat,Lon] = utm2deg(x,y,utmzone);
    
    xcomp=(Lat-NTR_Lat).^2;
    ycomp=(Lon-NTR_Lon).^2;
    d=sqrt(xcomp+ycomp);
    [~,I]=min(abs(d));
    
    NTR_name1=NTR_Lookup{I,3};

    if ii==1;
        load(['E:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\' NTR_name1 '.mat']);
        NTR_name2=NTR_name1;
    elseif ~strcmp(NTR_name2,NTR_name1);
        load(['E:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\' NTR_name1 '.mat']);
        NTR_name2=NTR_name1;
    else
    end
        
    
    NTR_index=NTR_Lookup{I,4};
    

    NTR_ts=NTR.WL(NTR_index,:);
    
    
    %Pick a fill_slope
    xcomp=(Lat-Fill_Slopes(:,1)).^2;
    ycomp=(Lon-Fill_Slopes(:,2)).^2;
    d=sqrt(xcomp+ycomp);
    [~,I]=min(abs(d));    
    fslope=Fill_Slopes(I,:);
    
    
   
    
    %Determine Slope and Runup formulation
    maxtide=max(tide)+max(NTR_ts);
     tide=tide+NTR_ts';
     clear NTR_ts
    [SlopeInfo]=SlopeMethod(wave_profile,test_con,Wavecon,cut,ESIID,Rey,maxtide,fslope);
    [RunupMethod]=Runup_Method_v1(SlopeInfo,wave_profile,Rey,fslope);
    
    StockTemp=RunupMethod.StockSlope;
    
    %% Now loop through the subprofiles and get their slopes as well
    %If the profle is 1, then do the first 5 in the folder
    %If at the end, do the last 5
    %Otherwise, look for the files that are 5 before and 5 after
    minordir=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_minor_profiles_wave_model_edited'];
    zip=[5:9,1:5];
    Rey=[];
    if ii==1;
        vsl=6:10;
    elseif ii==ii_grid_array(end);
        vsl=1:5;
    else
        vsl=1:10;
    end
        
    for jj=vsl;
        if jj<6;
            minorname=[minordir,filesep,'profile_', dec2base(ii-1,10,4),'_',num2str(zip(jj)),'.mat'];
        else
            minorname=[minordir,filesep,'profile_', dec2base(ii,10,4),'_',num2str(zip(jj)),'.mat'];
        end
        load(minorname);
        

        if wave_profile.msl_Lpos>300
            wave_profile.depth=wave_profile.depth(wave_profile.msl_Lpos-1:end);
            wave_profile.xpf=wave_profile.xpf(wave_profile.msl_Lpos-1:end);
            wave_profile.ypf=wave_profile.ypf(wave_profile.msl_Lpos-1:end);
            wave_profile.L=wave_profile.L(wave_profile.msl_Lpos-1:end);
            wave_profile.Lcoast= wave_profile.Lcoast-wave_profile.msl_Lpos-1;
            wave_profile.mhw_Lpos=wave_profile.mhw_Lpos-wave_profile.msl_Lpos+2;
            wave_profile.msl_Lpos=wave_profile.msl_Lpos-wave_profile.msl_Lpos+2;
        end      
        if isempty(wave_profile.depth);
            continue
        end
        [SlopeInfo]=SlopeMethod(wave_profile,test_con,Wavecon,cut,ESIID,Rey,maxtide,fslope);
        [RunupMethod]=Runup_Method_v1(SlopeInfo,wave_profile,Rey,fslope);    
        StockTemp=[StockTemp;RunupMethod.StockSlope];
              
    end
   
    
        StockSlope=mean(StockTemp);
    
    SlopeArray=[SlopeArray;StockSlope];
       

    
end
    name=[DirOut filesep NAME_GRD_transects '_StockSlopes']; 
    save(name,'SlopeArray'); 












    