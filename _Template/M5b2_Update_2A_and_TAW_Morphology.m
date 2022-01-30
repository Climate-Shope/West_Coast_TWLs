%For errors with 2A or rocky platform profiles, this script atempts to
%fix inconsistencies 

clear; close all; run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Douglas'; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI', filesep, NAME_GRD_transects, '_ShoreType.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_SeaIntersect.mat']);
load([directories.timeseries filesep 'ii_grid_array']);
if ~exist([directories.slope],'dir'), mkdir(directories.slope); end

utmzone='10 N';


%Load in the Precalculated regional Stockdon Slopes not included with these scripts
Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);


MorphDir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'Morphology_7_30_2019'];
MorphDir_3=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'Morphology_9_10_2019'];
if ~exist(MorphDir_3,'dir'), mkdir(MorphDir_3); end
cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');


%% First Flag the Profiles where there could be multiple profile types 
dirout      = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\Profile_2A']; if ~exist([dirout],'dir'), mkdir(dirout); end
file_esi    = ['F:\West_Coast_TWL_Hazards\01_Data\WashOre_ESI_Shapefiles' filesep 'Wash_Ore_ESI_UTM10_v2.shp'];
S=shaperead(file_esi);
load(['F:\West_Coast_TWL_Hazards\_STEP\_', NAME_GRD_transects,filesep, 'Transects2.mat']);

ESI_Array=[];
g=[];

for ii=1:length(Transects2.xpf);
     xpf=Transects2.xpf(ii,:);
     ypf=Transects2.ypf(ii,:);
     L1=[xpf;ypf];
     
     %if the profile read NaN, output NaNs
     if isnan(xpf(1));
           tempx=NaN;
           tempy=NaN;
           esi=NaN;
           Land=NaN;
           Sea=NaN;
           Trans=ii;
           
           g={Trans,tempx,tempy,esi,Land,Sea,NaN};
           ESI_Array=[ESI_Array;g];
     else
         
       for jj=1:numel(S);
       x=S(jj).X(1:end-1); 
       y=S(jj).Y(1:end-1);
       L2=[x;y];
       
       H1=mean(L1');
       H2=mean(L2');
       
       xcomp=(H1(1)-H2(1))^2;
       ycomp=(H1(2)-H2(2))^2;
       tdistance=sqrt(xcomp+ycomp);
       
       if tdistance > 5000
       continue
       else
       P=InterX(L1,L2);
       
       
       
       if isempty(P)
           continue
       else
           
           
           
           %Keep the intersection that is closest to the
           %offshore point
           xcomp=(P(1,:)-xpf(1)).^2;
           ycomp=(P(2,:)-ypf(1)).^2;
           distance=sqrt(xcomp+ycomp);
           [~,I]=min(distance);
           
            tempx=P(1,I);
           tempy=P(2,I);
           esi=S(jj).SEAWARD_SH;  
           k=strfind(esi,'2A');
           
            
           if isempty(k)
              continue 
           end
           Land=S(jj).LANDWARD_S;
           Sea=S(jj).SEAWARD_SH;
           Trans=ii;
           
            g={Trans,tempx,tempy,esi,Land,Sea,NaN};
           ESI_Array=[ESI_Array;g];
       end
       end
       end
     end
end


% Now from the ESI Array, flag the mixed Shorelines
ESI_2A=ESI_Array;

save([dirout filesep 'ESI_2A.mat'],'ESI_2A');


%%
% Load the potential 2A ESI data
load(['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\Profile_2A' filesep 'ESI_2A.mat']);
if ~isempty(ESI_2A) %flag if ESI2A is not empty
Flig=ESI_2A(:,1);
Flig=cell2mat(Flig);
else
  Flig=[];  
end

%load the Profile 
for ii=ii_grid_array;
        try
          load([MorphDir, filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat']);
        catch
          continue
        end
    
if strcmp(SlopeInfo.ESI,'2A');
    if ~isempty(SlopeInfo.berm_toe) && isempty(SlopeInfo.berm_width)
        SlopeInfo.berm_toe=[];
        
    elseif ~isempty(SlopeInfo.berm_width) && SlopeInfo.toe_loc<2 && SlopeInfo.toe_onshore<SlopeInfo.berm_loc+SlopeInfo.berm_width
        SlopeInfo.toe_loc=SlopeInfo.berm_loc+SlopeInfo.berm_width;
        SlopeInfo.toe_onshore=SlopeInfo.toe_loc;
        SlopeInfo.toe_ele=SlopeInfo.depth(SlopeInfo.toe_loc);
    elseif ~isempty(SlopeInfo.berm_width) && SlopeInfo.toe_loc<2 && SlopeInfo.toe_onshore>SlopeInfo.berm_loc+SlopeInfo.berm_width
        SlopeInfo.toe_loc=SlopeInfo.toe_onshore;
        SlopeInfo.toe_onshore=SlopeInfo.toe_loc;
        SlopeInfo.toe_ele=SlopeInfo.depth(SlopeInfo.toe_loc);
    end
end

if strcmp(SlopeInfo.Runup_Method.Method,'TAW') && SlopeInfo.toe_loc==1 || strcmp(SlopeInfo.Runup_Method.ExceedMethod,'TAW') && SlopeInfo.toe_loc==1 ;
    SlopeInfo.toe_ele=-1;
    SlopeInfo.toe_onshore=-1;
    SlopeInfo.toe_ele=SlopeInfo.MSL;
end


%% Now regardless, fix if the profile should be 2A and it isnt labeled as
% such 

if ismember(ii, Flig) && SlopeInfo.toe_loc>1;
    if strcmp(SlopeInfo.ESI,'3A') || strcmp(SlopeInfo.ESI,'4') || strcmp(SlopeInfo.ESI,'5') || strcmp(SlopeInfo.ESI,'6A')
    depth5=SlopeInfo.depth;
    result = DouglasPeucker([[1:length(depth5)];depth5],0.5);
    depth5=interp1(result(1,:),result(2,:),[1:length(depth5)]);
    depth5=depth5(1:SlopeInfo.toe_loc);
    tmp=detrend(depth5);
    [~,tmp]=max(tmp);
      bench_loc=tmp;
    bermwidth=SlopeInfo.toe_loc-bench_loc;
    bench_toe_loc=SlopeInfo.MHW_loc;
    if bench_toe_loc<0;
        bench_toe_loc=1;
    end
    
    if bench_loc==SlopeInfo.toe_loc;
        bench_loc=bench_toe_loc;
         bermwidth=SlopeInfo.toe_loc-bench_loc;
    end
 
    %%%Now add a check to ensure that the berm is roughly less than 1/15
    %%%slope as approximate by FEMA TAW guidelines
    if bench_loc>0;
        Clope=(SlopeInfo.depth(SlopeInfo.toe_loc)-SlopeInfo.depth(bench_loc))/(bermwidth);
    else
        Clope=(SlopeInfo.depth(SlopeInfo.toe_loc)-SlopeInfo.depth(1))/(SlopeInfo.toe_loc-1);
    end
    
    %if Clope>(1/15); %If the slope is not reasonaly horizontal, then wont consider it a berm condition roughly from the TAW guidelines
    if Clope>(1/12.5);
        bench_loc=[]; 
        bermwidth=[];
        
    end
    
    if SlopeInfo.toe_loc==1; %There cant really be a berm that can be measured in this case
        bench_loc=[]; 
        bermwidth=[]; 
        bench_toe_loc=[];
    end
    
    if ~isempty(bench_toe_loc) && isempty(bench_loc) && bench_toe_loc > 0;
        bench_loc=bench_toe_loc;
        bermwidth=SlopeInfo.toe_loc-bench_loc; 
        Clope=(SlopeInfo.depth(SlopeInfo.toe_loc)-SlopeInfo.depth(bench_loc))/(bermwidth);
    if Clope>(1/15) || bermwidth<2; %If the slope is not reasonaly horizontal, then wont consider it a berm condition roughly from the TAW guidelines
        bench_loc=[]; 
        bermwidth=[];
        
    end
    end
    
    if ~isempty(SlopeInfo.toe_loc)
    if bench_loc>=SlopeInfo.toe_loc
        bench_loc=[]; 
        bermwidth=[];
        bench_toe_loc=[];
    end
    end 
    
    SlopeInfo.berm_loc=bench_loc;
    SlopeInfo.berm_width=bermwidth;
    SlopeInfo.berm_toe=bench_toe_loc;
    
    if ~isempty(SlopeInfo.berm_width)
        gg=SlopeInfo.berm_loc-SlopeInfo.berm_toe;
        clope=(SlopeInfo.toe_ele-SlopeInfo.depth(SlopeInfo.berm_loc))/(SlopeInfo.toe_loc-SlopeInfo.berm_loc);
        ff=SlopeInfo.toe_loc-SlopeInfo.berm_loc;
        
        if gg>5 || abs(clope)>(1/12.5) || ff>25
           
            SlopeInfo.berm_loc=[]; 
            SlopeInfo.berm_width=[];
            SlopeInfo.berm_toe=[];            
            
        end
    end
    
    if ~isempty(SlopeInfo.berm_toe) && isempty(SlopeInfo.berm_width)
        SlopeInfo.berm_toe=[];
        
    elseif ~isempty(SlopeInfo.berm_width) && SlopeInfo.toe_loc<2 && SlopeInfo.toe_onshore<SlopeInfo.berm_loc+SlopeInfo.berm_width
        SlopeInfo.toe_loc=SlopeInfo.berm_loc+SlopeInfo.berm_width;
        SlopeInfo.toe_onshore=SlopeInfo.toe_loc;
        SlopeInfo.toe_ele=SlopeInfo.depth(SlopeInfo.toe_loc);
    elseif ~isempty(SlopeInfo.berm_width) && SlopeInfo.toe_loc<2 && SlopeInfo.toe_onshore>SlopeInfo.berm_loc+SlopeInfo.berm_width
        SlopeInfo.toe_loc=SlopeInfo.toe_onshore;
        SlopeInfo.toe_onshore=SlopeInfo.toe_loc;
        SlopeInfo.toe_ele=SlopeInfo.depth(SlopeInfo.toe_loc);
    end
    SlopeInfo.ESI='2A';
    end
    
end



    
 save([MorphDir_3, filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat'],'SlopeInfo');  
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
