%M7 Fix the extreme runup/overtopping values when TAW has calculated
%unrealistic runup values after initial cleanup attempt

clear all; close all; warning off
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m')  

%%
set(0, 'DefaultFigureColor','White','DefaultFigurePaperPositionMode','auto')
addpath(['F:\West_Coast_TWL_Hazards\_STEP\Extreme Value_Scripts\alt_stats'])

%% Initiate directories and laod in relevant data  
NAME_GRD_transects  = 'Douglas'; 

run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
addpath('F:\West_Coast_TWL_Hazards\_STEP\Extreme Value_Scripts\alt_stats\subCodes');
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI', filesep, NAME_GRD_transects, '_ShoreType.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_SeaIntersect.mat']);
load([directories.timeseries filesep 'ii_grid_array']);


utmzone='10 N';

Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);

MorphDir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'Morphology_9_10_2019'];
TWL_Dir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_9_16_2019'];
TWL_Dir_mid=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_update_10_12_2019'];
DirOut = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\extreme_value_analysis_update_10_14_2019'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end
DirOut_TWL = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_update_10_14_2019'];
if ~exist(DirOut_TWL,'dir'), mkdir(DirOut_TWL); end
return_period_array = [1.001;2;5;10;20;25;50;100;250;500];
load([directories.timeseries filesep 'time.mat'])
Count=0;

%% Determine a Cutoff TWL value
Socal=NAME_GRD_transects;
Val=[];
for kk=1;

load(['F:\West_Coast_TWL_Hazards\03_Results\Return_Levels_10_12_2019\' Socal '_RP10.mat']);
v=[RP10(:,3) RP10(:,16)];%Selecting RP10 as a general metric to help determien cutoff for abnormally high values
Val=[Val;v];
end
v=find(Val(:,2)==0);
Val(v,:)=[];
cutoff=mean(Val(:,1))+std(Val(:,1));%regional value for a threshold of high TAW TWLs above which need to investigate
clear RP10
load(['F:\West_Coast_TWL_Hazards\03_Results\Return_Levels_updated_10_12_2019\' NAME_GRD_transects '_RP10.mat']);

%% Evaluate individual profiles' extrema
for ii = ii_grid_array
    
    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);
    try
    %Load Extreme values 

        try
            AA=load(['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\extreme_value_analysis_update_10_12_2019\Extrema_GPD_',OBJECTID,'.mat']);
        catch
            AA=load(['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\extreme_value_analysis_update_10_12_2019\Extrema_prf',OBJECTID,'.mat']);
        end

       load([TWL_Dir_mid filesep 'TWL_Data_Transect_' num2str(ii) '.mat']); 
    SlopeInfo=TWL_Data.SlopeInfo;

    %Select correct toe
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


    if toe>1;
    Ej=SlopeInfo.depth(toe);
    else
    Ej=SlopeInfo.toe_ele; 
    end
    
    
%Flag TAW, only interested in these profiles as Stockdon does not return
%unrealistically large values 
    if isempty(TWL_Data.Composite_Slope);
      TAW_Flag=0;
   else
   if ~isempty(isnan(TWL_Data.Composite_Slope.Slope));
       TAW_Flag=1;
   else
       TAW_Flag=0;
   end
    end

    if TAW_Flag==0;
        continue
    end

    %See if R>3*Hmo at any point, if so, may be too large
    R=TWL_Data.Runup-TWL_Data.SWL;
    Hmo=(TWL_Data.DWL-Ej)*(0.78);
    Hmo(Hmo<1)=NaN;
    Too_Large=max(R./Hmo);
    if ~isnan(Too_Large)
        if Too_Large>3
            Too_Large=1;
        else
            Too_Large=0;
        end
    else
        Too_Large=0;
    end
    
    
    if Too_Large==1 || Ej<3 && AA.RLs(4)>cutoff; %if R>3*Hmo or the 10yr RP is greater than cutoff and the toe is < 3m

       Count=Count+1; 
       Flag(Count)=ii;
       zt(Count)=Ej;
    end
    catch
    end
    
end

%% Now with the flag, redo the TWL analysis, basically the same process as script M5c
cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');
%OKay first read the files in the directory 
listing = dir(['F:\West_Coast_TWL_Hazards\03_Results\Tides_Out\',NAME_GRD_transects]);
listing=size(listing);listing=listing(1)-2; 
tide_num=linspace(1,((listing-1)*10)+1,listing);

%Load in the NTR Lookup table
load('F:\West_Coast_TWL_Hazards\03_Results\SS and MMSL\NTR_Interp\NTR_Lookup');
%Grab the Output Coordinates
NTR_Lat=cell2mat(NTR_Lookup(:,1));
NTR_Lon=cell2mat(NTR_Lookup(:,2))*-1;


for ii=Flag;
    OBJECTID=dec2base(ii,10,4);
    
     try
    %Load wave timeseries and calculate the most energetic condition
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
    
    


    %% Pick the correct average Stockdon Slope for the profile, doesnt matter what the previous was determiend to be
    load([MorphDir filesep 'SlopeInfo_Transect_' num2str(ii) '.mat']);
    

    RunupMethod=SlopeInfo.Runup_Method;
    fslope=SlopeInfo.Runup_Method.StockSlope;
    fslope=[fslope,fslope,fslope];%Fill slope
    
   if SlopeInfo.toe_loc==1
      SlopeInfo.toe_loc=-1; 
   end
    
     if SlopeInfo.toe_loc>1;
    [Composite_Slope] = CompositeSlope(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
   
    else
    [Composite_Slope] = CompositeSlope_Plunging(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
     end
     
  %select correct toe      
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


if toe>1;
   Ej=SlopeInfo.depth(toe);
else
   Ej=SlopeInfo.toe_ele; 
end
 
 
%Same morphology changes as in M7 Tweak Runup version 3
    if Composite_Slope.wall_slope<36 && ~strcmp(RunupMethod.ESI,'6B') && strcmp(RunupMethod.Method,'Stockdon') %&& Ej<3
       RunupMethod.ExceedMethod='Stockdon';
    end
    
%% reassess morphology for the composite slope used for TAW for a few specific problem scenarios
    if Composite_Slope.wall_slope<45 && strcmp(RunupMethod.ESI,'6B') && SlopeInfo.depth(SlopeInfo.overtop_point)<13 %&& Ej<SlopeInfo.MSL
        SlopeInfo.overtop_point=SlopeInfo.maxi;
        zz=diff(diff(SlopeInfo.depth));
        zz(1:toe)=NaN;zz(SlopeInfo.overtop_point:end)=NaN; zz(SlopeInfo.depth>7)=NaN;
        [~,toe]=max(zz);
        SlopeInfo.toe_onshore=toe;
        SlopeInfo.toe_loc=toe;
        SlopeInfo.toe_ele=SlopeInfo.depth(toe);
        Ej=SlopeInfo.toe_ele;
        if SlopeInfo.maxi-SlopeInfo.CliffJun_Loc>5;
        zz=diff(diff(SlopeInfo.depth));
         zz(1:toe)=NaN;zz(SlopeInfo.overtop_point:end)=NaN; 
         [~,SlopeInfo.CliffJun_Loc]=min(zz);
         SlopeInfo.CliffJun_Ele=SlopeInfo.depth(SlopeInfo.CliffJun_Loc);
         SlopeInfo.CliffTop_Loc=SlopeInfo.CliffJun_Loc;
         SlopeInfo.CliffTop_Ele=SlopeInfo.CliffJun_Ele;
         SlopeInfo.overtop_point=SlopeInfo.CliffJun_Loc;
        end
     if SlopeInfo.toe_loc>1;
    [Composite_Slope] = CompositeSlope(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
      else
    [Composite_Slope] = CompositeSlope_Plunging(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
     end
    end
    
    
    if Composite_Slope.wall_slope<45 && SlopeInfo.depth(SlopeInfo.overtop_point)<13  && strcmp(RunupMethod.ExceedMethod,'TAW'); %&& Ej<3
       zz=diff(diff(SlopeInfo.depth));
        zz(SlopeInfo.overtop_point:end)=NaN; zz(SlopeInfo.depth>7)=NaN;
        [~,toe]=max(zz);
        SlopeInfo.toe_onshore=toe;
        SlopeInfo.toe_loc=toe;
        SlopeInfo.toe_ele=SlopeInfo.depth(toe);
        Ej=SlopeInfo.toe_ele;
        
         if SlopeInfo.maxi-SlopeInfo.CliffJun_Loc>5;
        zz=diff(diff(SlopeInfo.depth));
         zz(1:toe)=NaN;zz(SlopeInfo.overtop_point:end)=NaN; 
         [~,SlopeInfo.CliffJun_Loc]=min(zz);
         SlopeInfo.CliffJun_Ele=SlopeInfo.depth(SlopeInfo.CliffJun_Loc);
         SlopeInfo.CliffTop_Loc=SlopeInfo.CliffJun_Loc;
         SlopeInfo.CliffTop_Ele=SlopeInfo.CliffJun_Ele;
         end
    if SlopeInfo.toe_loc>1;
    [Composite_Slope] = CompositeSlope(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
     else
    [Composite_Slope] = CompositeSlope_Plunging(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);
     end
    end
    
    
    
    %now run TWLs through a more stringent TWL code where all Ib values > 8 are
    %selected for USACE method regardless of whther there is riprap
    %present or not
    [Runup,Runup_Mag,type,DWL] = CalcRunup_Fix6B_Ib8(Wavecon, SlopeInfo, RunupMethod, Composite_Slope,tide); 
 

 
    TWL_Data.SlopeInfo=SlopeInfo;
    TWL_Data.RunupMethod=RunupMethod;
    TWL_Data.Composite_Slope=Composite_Slope;
    TWL_Data.Runup=Runup; %<---- again this is TWL
    TWL_Data.Runup_Mag=Runup_Mag;
    TWL_Data.Depth=SlopeInfo.depth;
    TWL_Data.SWL=tide+SlopeInfo.MSL;
    TWL_Data.type=type;
    TWL_Data.DWL=DWL;
    name=[DirOut_TWL filesep 'TWL_Data_Transect_' num2str(ii)]; 
    save(name,'TWL_Data'); 
     catch
     end
end

clear TWL_Data Composite_Slope Wavecon SlopeInfo RunupMethod tide Ej NTR_ts

%% Now do an extreme value analysis same as version 3 approach
ii_ind=Flag;

Flag=[];
FlagDWL=[];
for ii = ii_ind
    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);
    fts=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_update_10_14_2019\TWL_Data_Transect_',num2str(ii),'.mat'];

    if ~exist(fts,'file'), continue, end
    TRANSID=ii;
    
    
    TWL = load(fts); 
    t=time;
    data=TWL.TWL_Data.Runup;
    %Distill to yearly maxima
    tvec=datevec(t);
    aa=datetime(tvec);
    TT = array2timetable(data,'RowTimes',aa);
    TT3 = retime(TT,'yearly','max');
    TT4=retime(TT,'monthly','max');
    Monthly_Data=TT4.data;
    Time_Monthly=datenum(TT4.Time);
    z_time=datevec(Time_Monthly);
    
    
    %Define Winter Here as Oct to Mar
    
    New_Vec=[];
    bb=find(z_time(:,1)==min(z_time(:,1)) & z_time(:,2)<4);
    vec=sortrows([Monthly_Data(bb) Time_Monthly(bb)],-1);
    vec=vec(1,:);
    New_Vec=[New_Vec;vec];
    kk=unique(z_time(:,1))';
    
    for ll=kk(2:end-1);
        bb=find(z_time(:,1)==ll-1 & z_time(:,2)>9);
        cc=find(z_time(:,1)==ll & z_time(:,2)<4);
        dd=[bb;cc];
        vec=sortrows([Monthly_Data(dd) Time_Monthly(dd)],-1);
        vec=vec(1,:);
        New_Vec=[New_Vec;vec];
    end
    
    bb=find(z_time(:,1)==max(kk) & z_time(:,2)>9);
    vec=sortrows([Monthly_Data(bb) Time_Monthly(bb)],-1);
    vec=vec(1,:);
    New_Vec=[New_Vec;vec];
    Yearly_Data=New_Vec(:,1);
    Time_Yearly=New_Vec(:,2);
    YEAR=kk';
    
    kk=ismember(data,Yearly_Data);
    type=TWL.TWL_Data.type(kk);
    
    
    clear TT TT3 tvec aa New_Vec vec
    try
    [res,GEV_RL,CDF_prob,CDF,GEV_conf]=GEV_generic_NO_FIG(Yearly_Data,YEAR,return_period_array);

    
    fout = [DirOut, filesep, 'Extrema_prf', dec2base(TRANSID,10,4),'.mat']; 
    m=matfile(sprintf('%s', fout),'writable',true);
    m.RP       = round(return_period_array); 
    m.RLs = GEV_RL; 
    m.FIT = 'GEV'; 
    m.CDF = CDF; 
    m.CDF_prob=CDF_prob;
    m.GEV_conf=GEV_conf;
    m.res=res;
    m.TRANSID = TRANSID;
    
    
    if GEV_conf(7)/GEV_RL(7)>0.35 || ~isreal(GEV_conf);
    Flag=[Flag;TRANSID];
    end
    
    %Now do the same for DWL 
    data=TWL.TWL_Data.DWL;
    tvec=datevec(t);
    aa=datetime(tvec);
    TT = array2timetable(data,'RowTimes',aa);
    TT3 = retime(TT,'yearly','max');
    
     TT4=retime(TT,'monthly','max');
    Monthly_Data=TT4.data;
    Time_Monthly=datenum(TT4.Time);
    z_time=datevec(Time_Monthly);
    
    %Define Winter Here as Oct to Mar
    
    New_Vec=[];
    bb=find(z_time(:,1)==min(z_time(:,1)) & z_time(:,2)<4);
    vec=sortrows([Monthly_Data(bb) Time_Monthly(bb)],-1);
    vec=vec(1,:);
    New_Vec=[New_Vec;vec];
    kk=unique(z_time(:,1))';
    
    for ll=kk(2:end-1);
        bb=find(z_time(:,1)==ll-1 & z_time(:,2)>9);
        cc=find(z_time(:,1)==ll & z_time(:,2)<4);
        dd=[bb;cc];
        vec=sortrows([Monthly_Data(dd) Time_Monthly(dd)],-1);
        vec=vec(1,:);
        New_Vec=[New_Vec;vec];
    end
    
    bb=find(z_time(:,1)==max(kk) & z_time(:,2)>9);
    vec=sortrows([Monthly_Data(bb) Time_Monthly(bb)],-1);
    vec=vec(1,:);
    New_Vec=[New_Vec;vec];
    Yearly_Data=New_Vec(:,1);
    Time_Yearly=New_Vec(:,2);
    YEAR=kk';

    
    kk=ismember(data,Yearly_Data);
    type=TWL.TWL_Data.type(kk);
    

    clear TT TT3 tvec aa  New_Vec vec

    [res_DWL,GEV_RL_DWL,CDF_prob_DWL,CDF_DWL,GEV_conf_DWL]=GEV_generic_NO_FIG(Yearly_Data,YEAR,return_period_array);
   
    m.RLs_DWL = GEV_RL_DWL; 
    m.CDF_DWL = CDF_DWL; 
    m.CDF_prob_DWL=CDF_prob_DWL;
    m.GEV_conf_DWL=GEV_conf_DWL;
    m.res_DWL=res_DWL;    
    
    if GEV_conf_DWL(7)/GEV_RL_DWL(7)>0.35 || ~isreal(GEV_conf_DWL);
    FlagDWL=[FlagDWL;TRANSID];
    end
    catch
      Flag=[Flag;TRANSID];
       FlagDWL=[FlagDWL;TRANSID];
    end
    
end


Investigate_Flag=[];

for ii=Flag';

    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);

    
    fts=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_update_10_14_2019\TWL_Data_Transect_',num2str(ii),'.mat'];

    if ~exist(fts,'file'), continue, end
    TRANSID=ii;
  
    TWL = load(fts); 
    t=time;
    data=TWL.TWL_Data.Runup;    
    tvec=datevec(t);
    aa=datetime(tvec);
    TT = array2timetable(data,'RowTimes',aa);
    TT2 = retime(TT,'daily','max');
    Daily_Data=TT2.data;
    Time_Data=datenum(TT2.Time);  
    
    kk=ismember(data,Daily_Data);
    type=TWL.TWL_Data.type(kk);
    
    tDclstr=3;
    [dataPks,tPks] = decluster([Time_Data Daily_Data],prctile(Daily_Data,90),tDclstr);
    
    kk=ismember(Daily_Data,dataPks);
    type=type(kk);
    
    thresh=prctile(Daily_Data,90);
    [GPD_RL,GPD_conf,CDF,CDF_prob,Num_Points,Percentile,res]=GPD_FIT_NO_FIG(dataPks,length(Daily_Data),return_period_array,type); 
    figure(1); 
    name=[DirOut filesep 'Profile_' num2str(TRANSID) '_GPD_Val_Plots.png'];
    saveas(gcf,name);
    figure(2); 
    name=[DirOut filesep 'Profile_' num2str(TRANSID) '_GPD_Return_Plots.png'];
    saveas(gcf,name);
    close all

    
    %Now decide which fit  to use 
    try
    First=load([DirOut, filesep, 'Extrema_prf', dec2base(TRANSID,10,4),'.mat']); 
    if First.GEV_conf(7)>GPD_conf(7) 
       fout = [DirOut, filesep, 'Extrema_GPD_prf', dec2base(TRANSID,10,4),'.mat']; 
       m=matfile(sprintf('%s', fout),'writable',true);
       m.RP  = round(return_period_array); 
       m.RLs = GPD_RL; 
       m.FIT = 'GPD'; 
       m.CDF = CDF; 
       m.CDF_prob=CDF_prob;
       m.GPD_conf=GPD_conf;
       m.res=res;
       m.TRANSID = TRANSID;
       m.Num_Points=Num_Points;
       m.Percentile=Percentile;
    end
    
    
    %Flag the Transects with bad fits regardless 
    if First.GEV_conf(7)>GPD_conf(7)
        if GPD_conf(7)/GPD_RL(7)>0.35;
           Investigate_Flag=[Investigate_Flag;TRANSID]; 
        end  
    else
        if First.GEV_conf(7)/First.RLs(7)>0.35;
           Investigate_Flag=[Investigate_Flag;TRANSID]; 
        end  
    end
    
    catch
       fout = [DirOut, filesep, 'Extrema_GPD_prf', dec2base(TRANSID,10,4),'.mat']; 
       m=matfile(sprintf('%s', fout),'writable',true);
       m.RP  = round(return_period_array); 
       m.RLs = GPD_RL; 
       m.FIT = 'GPD'; 
       m.CDF = CDF; 
       m.CDF_prob=CDF_prob;
       m.GPD_conf=GPD_conf;
       m.res=res;
       m.TRANSID = TRANSID;
       m.Num_Points=Num_Points;
       m.Percentile=Percentile;  
    end
    
end

for ii=FlagDWL';

    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);

    
    fts=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_update_10_14_2019\TWL_Data_Transect_',num2str(ii),'.mat'];

    if ~exist(fts,'file'), continue, end
    TRANSID=ii;
  
    TWL = load(fts); 
    t=time;
    data=TWL.TWL_Data.DWL;    
    tvec=datevec(t);
    aa=datetime(tvec);
    TT = array2timetable(data,'RowTimes',aa);
    TT2 = retime(TT,'daily','max');
    Daily_Data=TT2.data;
    Time_Data=datenum(TT2.Time);  
    
    kk=ismember(data,Daily_Data);
    type=TWL.TWL_Data.type(kk);
    
    tDclstr=3;
    [dataPks,tPks] = decluster([Time_Data Daily_Data],prctile(Daily_Data,90),tDclstr);
    
    kk=ismember(Daily_Data,dataPks);
    type=type(kk);
    
    thresh=prctile(Daily_Data,90);
    [GPD_RL,GPD_conf,CDF,CDF_prob,Num_Points,Percentile,res]=GPD_FIT_NO_FIG(dataPks,length(Daily_Data),return_period_array,type); 
    figure(1); 
    name=[DirOut filesep 'Profile_' num2str(TRANSID) '_GPD_Val_Plots.png'];
    saveas(gcf,name);
    figure(2); 
    name=[DirOut filesep 'Profile_' num2str(TRANSID) '_GPD_Return_Plots.png'];
    saveas(gcf,name);
    close all

    
    %Now decide which fit to use 
    try
    First=load([DirOut, filesep, 'Extrema_prf', dec2base(TRANSID,10,4),'.mat']); 
    if First.GEV_conf_DWL(7)>GPD_conf(7)
       fout = [DirOut, filesep, 'Extrema_GPD_prf_DWL_', dec2base(TRANSID,10,4),'.mat']; 
       m=matfile(sprintf('%s', fout),'writable',true);
       m.RP  = round(return_period_array); 
       m.RLs = GPD_RL; 
       m.FIT = 'GPD'; 
       m.CDF = CDF; 
       m.CDF_prob=CDF_prob;
       m.GPD_conf=GPD_conf;
       m.res=res;
       m.TRANSID = TRANSID;
       m.Num_Points=Num_Points;
       m.Percentile=Percentile;
    end
    
    
    %Flag the Transects with bad fits regardless 
    if First.GEV_conf_DWL(7)>GPD_conf(7)
        if GPD_conf_DWL(7)/GPD_RL(7)>0.35;
           Investigate_Flag_DWL=[Investigate_Flag;TRANSID]; 
        end  
    else
        if First.GEV_conf_DWL(7)/First.RLs(7)>0.35;
           Investigate_Flag_DWL=[Investigate_Flag;TRANSID]; 
        end  
    end
    catch
        fout = [DirOut, filesep, 'Extrema_GPD_prf_DWL_', dec2base(TRANSID,10,4),'.mat']; 
       m=matfile(sprintf('%s', fout),'writable',true);
       m.RP  = round(return_period_array); 
       m.RLs = GPD_RL; 
       m.FIT = 'GPD'; 
       m.CDF = CDF; 
       m.CDF_prob=CDF_prob;
       m.GPD_conf=GPD_conf;
       m.res=res;
       m.TRANSID = TRANSID;
       m.Num_Points=Num_Points;
       m.Percentile=Percentile; 
    end
    
end















