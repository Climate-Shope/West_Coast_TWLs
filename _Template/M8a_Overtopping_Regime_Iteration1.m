%Code to calculate the storm impact regime probability at each profile
%Preliminary Overtopping Regime
clear; close all; run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Douglas'; 
run('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');


%OUT
DirOut=['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Overtopping_10_2_2019'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end

%IN
TWLDir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'TWL_Output_9_16_2019'];
ExtDir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects filesep 'extreme_value_analysis_10_2_2019'];
load([directories.timeseries filesep 'ii_grid_array']);
load([directories.timeseries filesep 'time.mat'])
%4 categories:
%1 no impact on structure
%2 Collision
%3 Overtopping
%4 Inundation 


for ii=ii_grid_array;
    OBJECTID=dec2base(ii,10,4);
    try
    
   %First Load in the profile 
   load([TWLDir filesep 'TWL_Data_Transect_' num2str(ii) '.mat']); 
   %Now load in the Extrema
   name=[ExtDir filesep 'Extrema_GPD_prf' OBJECTID '.mat'];
   if exist(name,'file');
       load(name);
   else
     name=[ExtDir filesep 'Extrema_prf' OBJECTID '.mat']; 
     load(name);
   end
   name=[ExtDir filesep 'Extrema_GPD_prf_DWL_' OBJECTID '.mat'];
    if exist(name,'file');
       DWL=load(name);
    else
       DWL.RP  = RP; 
       DWL.RLs = RLs_DWL; 
       DWL.FIT = 'GEV'; 
       DWL.CDF = CDF_DWL; 
       DWL.CDF_prob=CDF_prob_DWL;
       DWL.GPD_conf=GEV_conf_DWL;
       DWL.res=res_DWL;
       DWL.TRANSID = TRANSID;
    end
   
   CDF_prob=CDF_prob(2:end-1);
   
   SlopeInfo=TWL_Data.SlopeInfo;
   TWL_ts=TWL_Data.Runup;
   DWL_ts=TWL_Data.DWL;
   clear TWL_Data
   MAT=[];

%Select appropriate toe elevation  
  
if SlopeInfo.toe_loc>=SlopeInfo.toe_onshore && SlopeInfo.toe_loc<SlopeInfo.overtop_point ;
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
  
  
  

  if toe<2
      toe_ele=SlopeInfo.toe_ele;
  else
      
  toe_ele=SlopeInfo.depth(toe);
  end
  
  
  yy=[toe_ele toe_ele];
  xx=[CDF_prob(1) CDF_prob(end)];
  
 %loop through return periods 
  for jj=1:length(RP);
      cdf1=CDF(jj,2:end-1);
      cdf2=DWL.CDF(jj,2:end-1);
      % no impact
      
      No_impact=find(cdf1<toe_ele);
      if isempty(No_impact)
         No_impact=0; 
      else
          
          [x,y]=curveintersect(xx,yy,CDF_prob,cdf1);
     
         if isempty(x) && cdf1(1)>toe_ele;
            No_impact=0;
         elseif isempty(x) && cdf1(end)<toe_ele;
             No_impact=1;
         else
            No_impact=x; 
         end

      end
      
      %Collision
      Collision1=find(cdf1>toe_ele);
      Collision2=find(cdf1<SlopeInfo.depth(SlopeInfo.overtop_point));
      
      if isempty(Collision1) %means we are in a full no impact regime
          Collision=0;
          
      else
          [x,y]=curveintersect(xx,yy,CDF_prob,cdf1);
          if isempty(x)
              Collision1=0;
          else
          Collision1=x;
          end
          
 
          [x,y]=curveintersect(xx,[SlopeInfo.depth(SlopeInfo.overtop_point) SlopeInfo.depth(SlopeInfo.overtop_point)],CDF_prob,cdf1);
          
          if isempty(x) && Collision1==0 && isempty(Collision2);
              Collision2=0;
          elseif isempty(x) && Collision1==0 && ~isempty(Collision2) && cdf1(end)<SlopeInfo.depth(SlopeInfo.overtop_point)
              Collision2=1;
          else
          
          Collision2=x;
          end
          
          if ~isempty(Collision2);
          Collision=Collision2-Collision1;
          else
            Collision=1-Collision1;  
          end
      end
      
      %Overtopping
      Overtop=find(cdf1>SlopeInfo.depth(SlopeInfo.overtop_point));
      if isempty(Overtop) 
          Overtop=0;   
      elseif length(Overtop)==length(cdf1)
          Overtop=1;
      else

          [x,y]=curveintersect(xx,[SlopeInfo.depth(SlopeInfo.overtop_point) SlopeInfo.depth(SlopeInfo.overtop_point)],CDF_prob,cdf1);
         Overtop=x;
         Overtop=1-Overtop;
      end
      
      %Inundation
      Inundation=find(cdf2>SlopeInfo.depth(SlopeInfo.overtop_point));
      if isempty(Inundation) 
          Inundation=0;
      else

          [x,y]=curveintersect(xx,[SlopeInfo.depth(SlopeInfo.overtop_point) SlopeInfo.depth(SlopeInfo.overtop_point)],CDF_prob,cdf2);
         if ~isempty(x)
          Inundation=x;
          Inundation=1-Inundation;
         else
           if cdf2(1)>SlopeInfo.depth(SlopeInfo.overtop_point)
               Inundation=1;
           else cdf2(end)<SlopeInfo.depth(SlopeInfo.overtop_point)
               Inundation=0;
               
           end
         end
          
          
      end
      
      %Now to get the appropriate Overtopping probability, 
      %subtract the inundation probability from the overtopping
      %probability
      Overtop=Overtop-Inundation;
      
      Probs=[No_impact,Collision,Overtop,Inundation];
      MAT=[MAT;Probs];
  end
  
 Impact.Probs=MAT; %10x4 matrix of impact probabilities
 Impact.Regime={'No_impact','Collision','Overtop','Inundation'};
 Impact.Regime_No=[1:4];
 Impact.RP=RP;
 Impact.SlopeInfo=SlopeInfo;
 Impact.OBJECTID=OBJECTID;
   

  
  
  
  %Now do DPY Analysis
  tvec=datevec(time);
  aa=datetime(tvec);
  TT = array2timetable(TWL_ts,'RowTimes',aa);
  TT2 = retime(TT,'daily','max');
  Daily_TWL=TT2.TWL_ts;
  TT = array2timetable(DWL_ts,'RowTimes',aa);
  TT2 = retime(TT,'daily','max');
  Daily_DWL=TT2.DWL_ts;
  
  
  N=length(find(Daily_TWL<toe_ele));
  C=length(find(find(Daily_TWL>toe_ele & Daily_TWL<SlopeInfo.depth(SlopeInfo.overtop_point))));
  O=length(find(Daily_TWL>SlopeInfo.depth(SlopeInfo.overtop_point)));
  I=length(find(Daily_DWL>SlopeInfo.depth(SlopeInfo.overtop_point)));
  O=O-I;
  
  Total=N+C+O+I;
  percents=[N/Total,C/Total,O/Total,I/Total];
  DPY=percents*365.25;
  
  Impact.DPY=DPY; %10x4 matrix of the number of days a profile experiences each regime in a given year
  Impact.DPY_Percents=percents;
  
   name=[DirOut filesep 'Impact_Transect_' num2str(ii) '.mat']; 
  save(name,'Impact');
  clear Impact MAT DPY DPY_Percents
    catch
 Impact.Probs=NaN;
 Impact.Regime={'No_impact','Collision','Overtop','Inundation'};
 Impact.Regime_No=[1:4];
 Impact.RP=RP;
 Impact.SlopeInfo=NaN;
 Impact.OBJECTID=OBJECTID;
 Impact.DPY=NaN;
  Impact.DPY_Percents=NaN;  
  name=[DirOut filesep 'Impact_Transect_' num2str(ii) '.mat']; 
  save(name,'Impact');   
        
    end
    
end