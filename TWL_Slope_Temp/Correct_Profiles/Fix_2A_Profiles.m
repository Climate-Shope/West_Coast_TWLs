%Fix values for 2A profiles, which are rocky bench/platform profiles which
% have specific morphology considerations not captured in base profile
% analysis code 

%% Define parameters
clear; close all; run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Santa_Cruz'; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI', filesep, NAME_GRD_transects, '_ShoreType.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_SeaIntersect.mat']);
load([directories.timeseries filesep 'ii_grid_array']);
if ~exist([directories.slope],'dir'), mkdir(directories.slope); end

%%  first read the files in the directory 
listing = dir(['F:\West_Coast_TWL_Hazards\03_Results\Tides_Out\',NAME_GRD_transects]);
listing=size(listing);listing=listing(1)-2; 
tide_num=linspace(1,((listing-1)*10)+1,listing);


%% Load in the NTR Lookup table
load('F:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\NTR_Lookup');
%Grab the Output Coordinates
NTR_Lat=cell2mat(NTR_Lookup(:,1));
NTR_Lon=cell2mat(NTR_Lookup(:,2))*-1;
utmzone='10 N';

%% Load in the Precalculated Beach Slopes
Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);
%%  Create the out directory 
DirOut=['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\TWL_Output_8_7_2019'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end
%% For testing purposes, create a random array of transects to run
tgg=randi(max(ii_grid_array),[9,1]);
tgg=[1;sort(tgg)]';

MorphDir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'Morphology_7_30_2019'];
cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');

%load the Profile 
for ii=ii_grid_array;
    try
     load([MorphDir, filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat']);
    catch
        continue
    end
    
%% Only care about 2A profiles in this script
if strcmp(SlopeInfo.ESI,'2A');
    if ~isempty(SlopeInfo.berm_width)
        gg=SlopeInfo.berm_loc-SlopeInfo.berm_toe;
        clope=(SlopeInfo.toe_ele-SlopeInfo.depth(SlopeInfo.berm_loc))/(SlopeInfo.toe_loc-SlopeInfo.berm_loc);
        ff=SlopeInfo.toe_loc-SlopeInfo.berm_loc;
        
        %Basically if these limits are exceeded, likely not a true Berm
        %shape so remove relevant parameters 
        if gg>5 || abs(clope)>(1/12.5) || ff>25
           
            SlopeInfo.berm_loc=[]; 
            SlopeInfo.berm_width=[];
            SlopeInfo.berm_toe=[];
        
            save([MorphDir, filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat'],'SlopeInfo');
            
        end
    end

            
          %Once Saved, time to Rerun the TWL Calc
          OBJECTID=dec2base(ii,10,4);
          Wavecon=load([directories.timeseries filesep 'WD4R_prf' dec2base(ii,10,4) '.mat']);
       
    
    if length(Wavecon.z0)>1;
        Wavecon.z0=Wavecon.z0(ii);
        Wavecon.x0=Wavecon.x0(ii);
        Wavecon.y0=Wavecon.y0(ii);
    end
    
    Wavecon.T(Wavecon.T<0)=0.1;
    
    [EI,idx]=max(Wavecon.Hs);
    test_con.Hs=Wavecon.Hs(idx);test_con.T=Wavecon.T(idx);
    
    wdepth=Wavecon.z0(1);
    if wdepth>10;
    [Ho,HoTwo]=BackCalcHo(test_con.T,wdepth,test_con.Hs);
    test_con.Hs=Ho;
    end
    
    %now load the profile 
    
    load([MorphDir, filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat']);
    ESIID=SlopeInfo.ESI;
    Rey=SlopeInfo.Rey;


    %Calc Tides
    %Find closest tidal output
    [~,temp]=min(abs(tide_num-ii));
    temp=tide_num(temp(1));
    load(['F:\West_Coast_TWL_Hazards\03_Results\Tides_Out\',NAME_GRD_transects,filesep, 'Tide_Series_',num2str(temp)]);
    tide=Tide_Series.WL;
    clear Tide_Series
    
    %Find and Load the appropriate NTR Dataset
    
    %Load the first point in the profile 
    x=SlopeInfo.xpf(1);
    y=SlopeInfo.ypf(1);
    
    [Lat,Lon] = utm2deg(x,y,utmzone);
    
    xcomp=(Lat-NTR_Lat).^2;
    ycomp=(Lon-NTR_Lon).^2;
    d=sqrt(xcomp+ycomp);
    [~,I]=min(abs(d));
    
    NTR_name1=NTR_Lookup{I,3};

    if ii==1;
        load(['F:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\' NTR_name1 '.mat']);
        NTR_name2=NTR_name1;
        
    elseif ~exist('NTR_name2','var')
        load(['F:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\' NTR_name1 '.mat']);
        NTR_name2=NTR_name1;
        
        
    elseif ~strcmp(NTR_name2,NTR_name1);
        load(['F:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\' NTR_name1 '.mat']);
        NTR_name2=NTR_name1;
    else
    end
        
    
    NTR_index=NTR_Lookup{I,4};
    
    NTR_ts=NTR.WL(NTR_index,:);
        tide=tide+NTR_ts';
    
    


    %% Pick the correct average Beach Slope for the profile, doesnt matter what the previous was determined to be
    load([MorphDir filesep 'SlopeInfo_Transect_' num2str(ii) '.mat']);

    RunupMethod=SlopeInfo.Runup_Method;
    fslope=SlopeInfo.Runup_Method.StockSlope;
    fslope=[fslope,fslope,fslope];
    
    
     if SlopeInfo.toe_loc>0;
        [Composite_Slope] = CompositeSlope(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
     else
        [Composite_Slope] = CompositeSlope_Plunging(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
     end
    
    [Runup,Runup_Mag,type] = CalcRunup_v7(Wavecon, SlopeInfo, RunupMethod, Composite_Slope,tide);
    TWL_Data.SlopeInfo=SlopeInfo;
    TWL_Data.RunupMethod=RunupMethod;
    TWL_Data.Composite_Slope=Composite_Slope;
    TWL_Data.Runup=Runup;
    TWL_Data.Runup_Mag=Runup_Mag;
    TWL_Data.Depth=SlopeInfo.depth;
    TWL_Data.SWL=tide+SlopeInfo.MSL;
    TWL_Data.type=type;
    name=[DirOut filesep 'TWL_Data_Transect_' num2str(ii)]; 
    save(name,'TWL_Data');  
end
    
    
    end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
