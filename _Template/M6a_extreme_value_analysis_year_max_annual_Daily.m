%Calculate Extreme value retun periods for TWL and DWL data

clear all; close all; warning off
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m')  

%%
set(0, 'DefaultFigureColor','White','DefaultFigurePaperPositionMode','auto')
addpath(['F:\West_Coast_TWL_Hazards\_STEP\Extreme Value_Scripts\alt_stats'])

%% CONFIG 
NAME_GRD_transects  = 'Douglas'; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
addpath('F:\West_Coast_TWL_Hazards\_STEP\Extreme Value_Scripts\alt_stats\subCodes');
%% LOAD TRANSECTS 
load(['F:\West_Coast_TWL_Hazards\_STEP\_', NAME_GRD_transects,filesep, 'Transects2.mat']); % loads transects and Ntr
load([directories.timeseries filesep 'ii_grid_array']);
Ntr = numel(Transects2.xpf(:,1));
%% Out
dirout = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\extreme_value_analysis_10_2_2019']; if ~exist(dirout,'dir'), mkdir(dirout); end

%%
VISIBILITY = 'on'; 
OVERWRITE = 1; 
MAKE_PLOT = 1; 

%%
% EXTREME ANALYSIS PARAMETERS 

return_period_array = [1.001;2;5;10;20;25;50;100;250;500]; %return period method cannot output a simple annual result, but 1.001 can be computed, whcih is close enough

% load time 
load([directories.timeseries filesep 'time.mat'])

% -------------------------------------------------------------------------
RUN_INDS = [ii_grid_array] ; % SUBSET TO RUN 
num = sort(randi(ii_grid_array(end),500,1)); %subset to generate diagnostic figures
c=find(ismember(num,ii_grid_array));
Fig_INDS=num(c)';
Fig_INDS=0;%remove line to generate figures
Flag=[];
FlagDWL=[];

for ii = ii_grid_array

    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);


    %% load time series of TWL
    
    fts=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_9_16_2019\TWL_Data_Transect_',num2str(ii),'.mat']; %Modify to appropriate directory of TWL output files

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
        
    %Fist try annual maxima GEV method
    if ismember(ii,Fig_INDS)
    [res,GEV_RL,CDF_prob,CDF,GEV_conf]=GEV_generic_FIG(Yearly_Data,YEAR,return_period_array,type);
    
    figure(1); 
    name=[dirout filesep 'Profile_' num2str(TRANSID) '_GEV_Val_Plots.png'];
    saveas(gcf,name);
    figure(2); 
    name=[dirout filesep 'Profile_' num2str(TRANSID) '_GEV_Return_Plots.png'];
    saveas(gcf,name);
    close all
    else
    [res,GEV_RL,CDF_prob,CDF,GEV_conf]=GEV_generic_NO_FIG(Yearly_Data,YEAR,return_period_array); %If no figure output is desired    
    end
    
    fout = [dirout, filesep, 'Extrema_prf', dec2base(TRANSID,10,4),'.mat']; 
    m=matfile(sprintf('%s', fout),'writable',true);
    m.RP       = round(return_period_array); 
    m.RLs = GEV_RL; 
    m.FIT = 'GEV'; 
    m.CDF = CDF; 
    m.CDF_prob=CDF_prob;
    m.GEV_conf=GEV_conf;
    m.res=res;
    m.TRANSID = TRANSID;
    
    
    if GEV_conf(8)/GEV_RL(8)>0.35 || ~isreal(GEV_conf);
    Flag=[Flag;TRANSID]; %If listed in FLAG, profile will be reassessed using GPD method, threshold of 0.35 determiend by testing multiple profiles
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
   
    if GEV_conf_DWL(8)/GEV_RL_DWL(8)>0.35 || ~isreal(GEV_conf_DWL);
        FlagDWL=[FlagDWL;TRANSID];
    end
    catch
      Flag=[Flag;TRANSID];
      FlagDWL=[FlagDWL;TRANSID];%If listed in FLAGDWL, profile will be reassessed using GPD method
    end
    
    
end

%% Now to select those with Flag and see if GPD is a better fit

Investigate_Flag=[];

for ii=Flag';

    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);

    
    fts=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_9_16_2019\TWL_Data_Transect_',num2str(ii),'.mat'];

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
    [GPD_RL,GPD_conf,CDF,CDF_prob,Num_Points,Percentile,res]=GPD_FIT_FIG(dataPks,length(Daily_Data),return_period_array,type); 
    figure(1); 
    %Save output figures for diagnostics
    name=[dirout filesep 'Profile_' num2str(TRANSID) '_GPD_Val_Plots.png'];
    saveas(gcf,name);
    figure(2); 
    name=[dirout filesep 'Profile_' num2str(TRANSID) '_GPD_Return_Plots.png'];
    saveas(gcf,name);
    close all

    
    %Now decide which fit to use 
    try
    First=load([dirout, filesep, 'Extrema_prf', dec2base(TRANSID,10,4),'.mat']); 
    if First.GEV_conf(8)>GPD_conf(8) % If GPD has a smaller CI, save it, otherwise use GEV
       fout = [dirout, filesep, 'Extrema_GPD_prf', dec2base(TRANSID,10,4),'.mat']; 
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
    if First.GEV_conf(8)>GPD_conf(8)
        if GPD_conf(8)/GPD_RL(8)>0.35;
           Investigate_Flag=[Investigate_Flag;TRANSID]; 
        end  
    else
        if First.GEV_conf(8)/First.RLs(8)>0.35;
           Investigate_Flag=[Investigate_Flag;TRANSID]; 
        end  
    end
    
    catch
       fout = [dirout, filesep, 'Extrema_GPD_prf', dec2base(TRANSID,10,4),'.mat']; 
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
save([dirout, filesep, 'Investigate_Flag.mat'],'Investigate_Flag');
save([dirout, filesep, 'Flag.mat'],'Flag');



Investigate_Flag=[];
Investigate_Flag_DWL=[];

for ii=FlagDWL';

    TRANSID = ii; 
    OBJECTID=dec2base(ii,10,4);

    
    fts=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\TWL_Output_9_16_2019\TWL_Data_Transect_',num2str(ii),'.mat'];

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
    [GPD_RL,GPD_conf,CDF,CDF_prob,Num_Points,Percentile,res]=GPD_FIT_FIG(dataPks,length(Daily_Data),return_period_array,type); 
    figure(1); 
    %Output figures for diagnostics
    name=[dirout filesep 'Profile_' num2str(TRANSID) '_GPD_Val_Plots.png'];
    saveas(gcf,name);
    figure(2); 
    name=[dirout filesep 'Profile_' num2str(TRANSID) '_GPD_Return_Plots.png'];
    saveas(gcf,name);
    close all

    
    %Now decide which fit to use 
    try
    First=load([dirout, filesep, 'Extrema_prf', dec2base(TRANSID,10,4),'.mat']); 
    if First.GEV_conf_DWL(8)>GPD_conf(8)% If GPD has a namller CI, save it, otherwise use GEV
       fout = [dirout, filesep, 'Extrema_GPD_prf_DWL_', dec2base(TRANSID,10,4),'.mat']; 
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
    if First.GEV_conf_DWL(8)>GPD_conf(8)
        if GPD_conf_DWL(8)/GPD_RL(8)>0.35;
           Investigate_Flag_DWL=[Investigate_Flag;TRANSID]; 
        end  
    else
        if First.GEV_conf_DWL(8)/First.RLs(8)>0.35;
           Investigate_Flag_DWL=[Investigate_Flag;TRANSID]; 
        end  
    end
    catch
        fout = [dirout, filesep, 'Extrema_GPD_prf_DWL_', dec2base(TRANSID,10,4),'.mat']; 
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
save([dirout, filesep, 'Investigate_Flag_DWL.mat'],'Investigate_Flag_DWL');
save([dirout, filesep, 'FlagDWL.mat'],'FlagDWL');


