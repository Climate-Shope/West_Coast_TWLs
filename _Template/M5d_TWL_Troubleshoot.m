%Troubleshoot some TWL profiles that did not yield results 
clear; close all; run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Douglas'; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI', filesep, NAME_GRD_transects, '_ShoreType.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_SeaIntersect.mat']);
load([directories.timeseries filesep 'ii_grid_array']);
if ~exist([directories.slope],'dir'), mkdir(directories.slope); end

%OKay first read the files in the directory 
listing = dir(['F:\West_Coast_TWL_Hazards\03_Results\Tides_Out\',NAME_GRD_transects]);
listing=size(listing);listing=listing(1)-2; 
tide_num=linspace(1,((listing-1)*10)+1,listing);


%Load in the NTR Lookup table
load('F:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\NTR_Lookup');
%Grab the Output Coordinates
NTR_Lat=cell2mat(NTR_Lookup(:,1));
NTR_Lon=cell2mat(NTR_Lookup(:,2))*-1;

utmzone='10 N';


%Load in the Precalculated Stockdon Slopes
Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);
% Create the out directory 
DirOut=['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\TWL_Output_9_16_2019'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end


MorphDir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'Morphology_9_10_2019'];
cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');

load([DirOut filesep 'Run_Success.mat'])
tgg=find(strcmp(Run_Success.Complete,'Fail')); %only look at profiles that failed to run
if isempty(tgg)
    return
end
for ii=tgg;

    OBJECTID=dec2base(ii,10,4);
    try
    %Load wave timeseries and calculate the most energetic condtion
    Wavecon=load([directories.timeseries filesep 'WD4R_prf' dec2base(ii,10,4) '.mat']);
       
    
    if length(Wavecon.z0)>1;
        Wavecon.z0=Wavecon.z0(ii);
        Wavecon.x0=Wavecon.x0(ii);
        Wavecon.y0=Wavecon.y0(ii);
    end
    
    Wavecon.T(Wavecon.T<0)=0.1;
    Wavecon.T(isnan(Wavecon.T))=0.1;
    Wavecon.Hs(isnan(Wavecon.Hs))=0.1;
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
    %Find closest
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
    
    


    %% Pick the correct precomputed average Stockdon Slope for the profile, doesnt matter what the previous was determiend to be
    load([MorphDir filesep 'SlopeInfo_Transect_' num2str(ii) '.mat']);
    

    RunupMethod=SlopeInfo.Runup_Method;
    fslope=SlopeInfo.Runup_Method.StockSlope;
    fslope=[fslope,fslope,fslope];
    
    
   if SlopeInfo.toe_onshore==-1 && SlopeInfo.toe_loc~=-1;
       if SlopeInfo.toe_ele>5
           SlopeInfo.toe_loc=-1;
           SlopeInfo.toe_ele=SlopeInfo.MSL;
       else
           SlopeInfo.toe_onshore=SlopeInfo.toe_onshore;
       end
       
   end
    
    if isnan(SlopeInfo.Runup_Method.StockSlope);
        a=find(strcmp(AllStockdon.Counties,NAME_GRD_transects));
        b=find(strcmp(AllStockdon.Names,['profile_' OBJECTID]));
        c=find(ismember(b,a));
        c=b(c);
        d=find(~isnan(AllStockdon.Slope_MeanMajors_Limited));
        [~,e]=min(abs(d-c));
        f=d(e);
        SlopeInfo.Runup_Method.StockSlope=AllStockdon.Slope_MeanMajors_Limited(f);
        fslope=SlopeInfo.Runup_Method.StockSlope;
        fslope=[fslope,fslope,fslope];
    end
    
    
    
    
 if SlopeInfo.toe_loc>=SlopeInfo.toe_onshore && SlopeInfo.toe_loc<SlopeInfo.overtop_point;
    toe=SlopeInfo.toe_loc;
    
    if toe>0;
    if SlopeInfo.depth(toe) > 6 && SlopeInfo.depth(SlopeInfo.toe_onshore)<6;
        toe=SlopeInfo.toe_onshore;
    end
    end
else
    toe=SlopeInfo.toe_onshore;
    
    if SlopeInfo.depth(toe) < 6 && SlopeInfo.depth(SlopeInfo.toe_onshore)>6;
        toe=SlopeInfo.toe_loc;
    end
end
    
    if toe>SlopeInfo.overtop_point
        zz=movingslope(SlopeInfo.depth(1:SlopeInfo.overtop_point),3);
        zz=movingslope(zz,3);
        [~,KK]=max(zz);
        if KK>3.5
            SlopeInfo.toe_loc=-1;
            SlopeInfo.toe_onshore=-1;
            SlopeInfo.toe_ele=SlopeInfo.MSL;
        else
           SlopeInfo.toe_loc=KK;
           SlopeInfo.toe_onshore=KK;
           SlopeInfo.toe_ele= SlopeInfo.depth(KK);
        end
    end
    
    
    
   if SlopeInfo.toe_loc==1
      SlopeInfo.toe_loc=-1; 
   end
   

    
     if SlopeInfo.toe_loc>1;
    [Composite_Slope] = CompositeSlope(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
   
    else
    [Composite_Slope] = CompositeSlope_Plunging(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
    end
    if ~isempty(Composite_Slope);
        pause_space=1;
    end
    
    [Runup,Runup_Mag,type,DWL] = CalcRunup_v7(Wavecon, SlopeInfo, RunupMethod, Composite_Slope,tide);

    
    
    TWL_Data.SlopeInfo=SlopeInfo;
    TWL_Data.RunupMethod=RunupMethod;
    TWL_Data.Composite_Slope=Composite_Slope;
    TWL_Data.Runup=Runup;
    TWL_Data.Runup_Mag=Runup_Mag;
    TWL_Data.Depth=SlopeInfo.depth;
    TWL_Data.SWL=tide+SlopeInfo.MSL;
    TWL_Data.type=type;
    TWL_Data.DWL=DWL;
    name=[DirOut filesep 'TWL_Data_Transect_' num2str(ii)]; 
    save(name,'TWL_Data');
    
    Fix_Success.OBJECTID{ii}=OBJECTID;
    Fix_Success.Complete{ii}='Pass';
    catch
        if exist([MorphDir filesep 'SlopeInfo_Transect_' num2str(ii) '.mat'],'file')
    Fix_Success.OBJECTID{ii}=OBJECTID;
    Fix_Success.Complete{ii}='Fail';    
        else
         Fix_Success.OBJECTID{ii}=OBJECTID;
        Fix_Success.Complete{ii}='Pass';       
        end
    end
    
    
end
    name=[DirOut filesep 'Fix_Success']; %For diagnostics
    save(name,'Fix_Success');
    