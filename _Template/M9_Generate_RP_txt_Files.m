%% Code to generate text files of Extreme value outputs, profile locations, water levels,and other critical information 

%Define Parameters
clear; close all; run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = {'San_Diego','Orange','Los_Angeles','Ventura','Santa_Barbara',...
    'Santa_Barbara_West','San_Luis_Obispo','Monterey','Santa_Cruz','San_Mateo','San_Francisco',...
    'Marin','Sonoma','Mendocino','Humboldt','Del_Norte','Curry','Coos','Douglas','Lane',...
    'Lincoln','Tillamook','Clatsop','Pacific','Grays_Harbor','Jefferson','Clallam'}; 
run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\run_init_directories.m') 
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
addpath('E:/West_Coast_TWL_Hazards/_STEP/_Compile');
UTM={'11 N','11 N','11 N','11 N','11 N','10 N','10 N','10 N','10 N',...
     '10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N',...
     '10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N','10 N'};
DirOut=['E:\West_Coast_TWL_Hazards\03_Results\Return_Levels_updated_03_16_2020'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end 
 
 for WestCoast=1:27;
  %DirTWL=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\TWL_Output_9_16_2019'];   
  DirOVR=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\Overtopping_10_2_2019'];   
  DirEXT=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\extreme_value_analysis_10_2_2019'];
  DirOVR_up=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\Overtopping_updated_10_14_2019'];   
  DirEXT_up=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\extreme_value_analysis_update_10_14_2019'];
  DirTWL=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\TWL_Output_9_16_2019'];
  DirTWL_up=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\TWL_Output_update_10_14_2019'];
  DirOVR_mid=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\Overtopping_updated_10_12_2019'];   
  DirEXT_mid=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\extreme_value_analysis_update_10_12_2019'];
  DirTWL_mid=['E:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\TWL_Output_update_10_12_2019'];
  load(['E:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\time_series_reconstructed',filesep,'ii_grid_array.mat']);
 
  MHW_Dir=['E:\West_Coast_TWL_Hazards\03_Results\MHW' filesep NAME_GRD_transects{WestCoast} '_MHW.mat'];
   
  MHW_file=load(MHW_Dir);
  
  
  Count=0;
  for ii= ii_grid_array;

    try
       if exist([DirEXT_up filesep 'Extrema_prf' dec2base(ii,10,4) '.mat'],'file')
               OVERTOP=load([DirOVR_up filesep 'Impact_Transect_' num2str(ii)]);
      
        try 
        Extrema_TWL=load([DirEXT_up filesep 'Extrema_GPD_prf' dec2base(ii,10,4)]);
        catch
        Extrema_TWL=load([DirEXT_up filesep 'Extrema_prf' dec2base(ii,10,4)]);    
        end
      
        try 
        Extrema_DWL=load([DirEXT_up filesep 'Extrema_GPD_prf_DWL_' dec2base(ii,10,4)]);
        catch
        Extrema_DWL=load([DirEXT_up filesep 'Extrema_prf' dec2base(ii,10,4)]);      
        end
       
       elseif exist([DirEXT_mid filesep 'Extrema_prf' dec2base(ii,10,4) '.mat'],'file')
               OVERTOP=load([DirOVR_mid filesep 'Impact_Transect_' num2str(ii)]);
      
        try 
        Extrema_TWL=load([DirEXT_mid filesep 'Extrema_GPD_prf' dec2base(ii,10,4)]);
        catch
        Extrema_TWL=load([DirEXT_mid filesep 'Extrema_prf' dec2base(ii,10,4)]);    
        end
      
        try 
        Extrema_DWL=load([DirEXT_mid filesep 'Extrema_GPD_prf_DWL_' dec2base(ii,10,4)]);
        catch
        Extrema_DWL=load([DirEXT_mid filesep 'Extrema_prf' dec2base(ii,10,4)]);      
        end 
           
       else
      %Load TWL  

      OVERTOP=load([DirOVR filesep 'Impact_Transect_' num2str(ii)]);
      
      try 
      Extrema_TWL=load([DirEXT filesep 'Extrema_GPD_prf' dec2base(ii,10,4)]);
      catch
      Extrema_TWL=load([DirEXT filesep 'Extrema_prf' dec2base(ii,10,4)]);    
      end
      
      try 
      Extrema_DWL=load([DirEXT filesep 'Extrema_GPD_prf_DWL_' dec2base(ii,10,4)]);
      catch
      Extrema_DWL=load([DirEXT filesep 'Extrema_prf' dec2base(ii,10,4)]);      
      end
      end
   %Get LAT LON using the MHW Contour
   U=OVERTOP.Impact.SlopeInfo.xpf(OVERTOP.Impact.SlopeInfo.MHW_loc);
   V=OVERTOP.Impact.SlopeInfo.ypf(OVERTOP.Impact.SlopeInfo.MHW_loc);
   
   [Lat,Lon] = utm2deg(U,V,UTM{WestCoast});

   
   try
       TWL_Data=load([DirTWL_up filesep 'TWL_Data_Transect_' num2str(ii)]);
   catch
       try
           TWL_Data=load([DirTWL_mid filesep 'TWL_Data_Transect_' num2str(ii)]);
       catch
       
        TWL_Data=load([DirTWL filesep 'TWL_Data_Transect_' num2str(ii)]);
       end
   end
   
   MHW=TWL_Data.TWL_Data.SlopeInfo.MHW;
   
   
   if isempty(TWL_Data.TWL_Data.Composite_Slope);
      TAW_Flag=0;
   else
   if ~isempty(isnan(TWL_Data.TWL_Data.Composite_Slope.Slope));
       TAW_Flag=1;
   else
       TAW_Flag=0;
   end
   end
       
       
   if Lon> -119 && WestCoast==5
       [Lat,Lon] = utm2deg(U,V,'10 N');
   end;
   
   if Lon> -5 
       continue
   end;
   
   if Extrema_TWL.RLs(4)>12
   Pause=1;
   end
    
   %Get_Toe
   SlopeInfo=OVERTOP.Impact.SlopeInfo;
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
    
    
   %Get xy of each WL point
   %TWL
   for gg=1:10;
       temp1=TWL_Data.TWL_Data.SlopeInfo.depth-Extrema_TWL.RLs(gg);
       temp2=find(temp1>0);
       if ~isempty(temp2)
           location=temp2(1);
       else
           location=length(TWL_Data.TWL_Data.SlopeInfo.depth);
       end
       
       TWL_Lat=TWL_Data.TWL_Data.SlopeInfo.ypf(location);
       TWL_Lon=TWL_Data.TWL_Data.SlopeInfo.xpf(location);
       [TWL_Lat_Deg,TWL_Lon_Deg] = utm2deg(TWL_Lon,TWL_Lat,UTM{WestCoast});
       TWL_Latitude(gg)=TWL_Lat_Deg;
       TWL_Longitude(gg)=TWL_Lon_Deg;
   end
  
   clear temp1 temp2 location TWL_Lat TWL_Lon TWL_Lat_Deg TWL_Lon_Deg
   
   %DWL
    for gg=1:10;
      if Extrema_DWL.RLs_DWL(gg)>Extrema_TWL.RLs(gg)
         Extrema_DWL.RLs_DWL(gg)=Extrema_TWL.RLs(gg);
      end
      
       temp1=TWL_Data.TWL_Data.SlopeInfo.depth-Extrema_DWL.RLs_DWL(gg);
       
       
       temp2=find(temp1>0);
       if ~isempty(temp2)
           location=temp2(1);
       else
           location=length(TWL_Data.TWL_Data.SlopeInfo.depth);
       end
       
       DWL_Lat=TWL_Data.TWL_Data.SlopeInfo.ypf(location);
       DWL_Lon=TWL_Data.TWL_Data.SlopeInfo.xpf(location);
       [DWL_Lat_Deg,DWL_Lon_Deg] = utm2deg(DWL_Lon,DWL_Lat,UTM{WestCoast});
       DWL_Latitude(gg)=DWL_Lat_Deg;
       DWL_Longitude(gg)=DWL_Lon_Deg;
   end
   
   clear temp1 temp2 location DWL_Lat DWL_Lon DWL_Lat_Deg DWL_Lon_Deg
   
   % get transect endpoint
   [EndLat,EndLon] = utm2deg(TWL_Data.TWL_Data.SlopeInfo.xpf(end),TWL_Data.TWL_Data.SlopeInfo.ypf(end),UTM{WestCoast});
    
   
   %Fix MHW Estimate and location
   
   
   if MHW>2.5
       %if greater than 2.5, need to find the closest point and use that
        X=TWL_Data.TWL_Data.SlopeInfo.xpf(1);
        Y=TWL_Data.TWL_Data.SlopeInfo.ypf(1);
        
        %Find distance
        xcomp=(MHW_file.MHW(:,1)-X).^2;
        ycomp=(MHW_file.MHW(:,2)-Y).^2;
        d=sqrt(xcomp+ycomp);
        [~,i]=min(d);
        MHW=MHW_file.MHW(1,3);
   end
   
   [MHW_Lat,MHW_Lon] = utm2deg(TWL_Data.TWL_Data.SlopeInfo.xpf(TWL_Data.TWL_Data.SlopeInfo.MHW_loc),TWL_Data.TWL_Data.SlopeInfo.ypf(TWL_Data.TWL_Data.SlopeInfo.MHW_loc),UTM{WestCoast});
   
   %toe and crest locations
   
   try
   [Toe_Lat,Toe_Lon] = utm2deg(TWL_Data.TWL_Data.SlopeInfo.xpf(TWL_Data.TWL_Data.SlopeInfo.toe_loc),TWL_Data.TWL_Data.SlopeInfo.ypf(TWL_Data.TWL_Data.SlopeInfo.toe_loc),UTM{WestCoast});
   catch
   Toe_Lat=MHW_Lat;
   Toe_Lon=MHW_Lon;
   end
   
   try
   [Cliff_Lat,Cliff_Lon] = utm2deg(TWL_Data.TWL_Data.SlopeInfo.xpf(TWL_Data.TWL_Data.SlopeInfo.overtop_point),TWL_Data.TWL_Data.SlopeInfo.ypf(TWL_Data.TWL_Data.SlopeInfo.overtop_point),UTM{WestCoast});
   catch 
   Cliff_Lat=EndLat;
   Cliff_Lon=EndLon;
   end
   
   
   %Fix all of the coordinates that are erroneous 
   %Lon
   a=find(Lon>-116);
   if ~isempty(a)
   b=Lon(a);
   c=Lat(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
   Lat(a)=y2;
   Lon(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end
   
   %EndLon
   a=find(EndLon>-116);
   if ~isempty(a)
   b=EndLon(a);
   c=EndLat(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
   EndLat(a)=y2;
   EndLon(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end
   
   %TWL_Longitude
   a=find(TWL_Longitude>-116);
   if ~isempty(a)
   
   b=TWL_Longitude(a);
   c=TWL_Latitude(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
     
  
   TWL_Latitude(a)=y2;
   TWL_Longitude(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end
 
   %DWL_Longitude
   a=find(DWL_Longitude>-116);
   if ~isempty(a)
   b=DWL_Longitude(a);
   c=DWL_Latitude(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
   DWL_Latitude(a)=y2;
   DWL_Longitude(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end 
 
   %Toe_Lon
   a=find(Toe_Lon>-116);
   if ~isempty(a)
   b=Toe_Lon(a);
   c=Toe_Lat(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
   Toe_Lat(a)=y2;
   Toe_Lon(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end 
   
   %Cliff_Lon
   a=find(Cliff_Lon>-116);
   if ~isempty(a)
   b=Cliff_Lon(a);
   c=Cliff_Lat(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
   Cliff_Lat(a)=y2;
   Cliff_Lon(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end 

   %MHW_Lon
   a=find(MHW_Lon>-116);
   if ~isempty(a)
   b=MHW_Lon(a);
   c=MHW_Lat(a);
   [x,y,blah] = deg2utm(c,b)
   y2=[];
   x2=[];
   for rr=1:length(x);
   [y1,x1] = utm2deg(x(rr),y(rr),'10 N')
   y2=[y2;y1];
   x2=[x2;x1];
   end
   MHW_Lat(a)=y2;
   MHW_Lon(a)=x2;
   clear a b c x y blah y2 x2 y1 x1
   end    
   
      
   %initialize each RP    
   for kk=1:10;
     if Extrema_DWL.RLs_DWL(kk)>Extrema_TWL.RLs(kk)
       Extrema_DWL.RLs_DWL(kk)=Extrema_TWL.RLs(kk);
     end
     
     
       tmp(kk,:)=[round(Lat,6),round(Lon,6),round(EndLat,6),round(EndLon,6),round(Extrema_TWL.RLs(kk),2),round(TWL_Latitude(kk),6),round(TWL_Longitude(kk),6),round(Extrema_DWL.RLs_DWL(kk),2),round(DWL_Latitude(kk),6),round(DWL_Longitude(kk),6),round(toe_ele,2),round(Toe_Lat,6), round(Toe_Lon,6),round(SlopeInfo.depth(SlopeInfo.overtop_point),2),round(Cliff_Lat,6),round(Cliff_Lon,6),MHW,round(MHW_Lat,6), round(MHW_Lon,6), OVERTOP.Impact.Probs(kk,:),round(OVERTOP.Impact.DPY,2)];

   end
   
   Count=Count+1; 
   %Separate for Appropriate RP 
   RP1(Count,:)=tmp(1,:);
   RP2(Count,:)=tmp(1,:);
   RP5(Count,:)=tmp(2,:);
   RP10(Count,:)=tmp(3,:);
   RP20(Count,:)=tmp(4,:);
   RP25(Count,:)=tmp(5,:);
   RP50(Count,:)=tmp(6,:);
   RP100(Count,:)=tmp(7,:);
   RP250(Count,:)=tmp(8,:);
   RP500(Count,:)=tmp(9,:);
   catch

   end
   
      
      
  end
  
  if WestCoast>16
     RP1=flipud(RP1);
     RP2=flipud(RP2); 
     RP5=flipud(RP5); 
     RP10=flipud(RP10); 
     RP20=flipud(RP20); 
     RP25=flipud(RP25); 
     RP50=flipud(RP50); 
     RP100=flipud(RP100); 
     RP250=flipud(RP250); 
     RP500=flipud(RP500); 
  end
  
  
  
  
  
  %Save the Matricies Separately
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP1.mat'];
  save(Name,'RP1');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP2.mat'];
  save(Name,'RP2');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP5.mat'];
  save(Name,'RP5');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP10.mat'];
  save(Name,'RP10');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP20.mat'];
  save(Name,'RP20');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP25.mat'];
  save(Name,'RP25');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP50.mat'];
  save(Name,'RP50');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP100.mat'];
  save(Name,'RP100');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP250.mat'];
  save(Name,'RP250');
  Name=[DirOut filesep NAME_GRD_transects{WestCoast} '_RP500.mat'];
  save(Name,'RP500');
  
  
  clear RP1 RP2 RP5 RP10 RP20 RP25 RP50 RP100 RP250 RP500 OVERTOP Extrema_TWL Extrema_DWL TWL_Data MHW TAW_Flag 
      
     
 end
 
 
 
 %% Now Compile all into one WestCoast Mat File
 RP1_Westcoast=[];
 RP2_Westcoast=[];
 RP5_Westcoast=[];
 RP10_Westcoast=[];
 RP20_Westcoast=[];
 RP25_Westcoast=[];
 RP50_Westcoast=[];
 RP100_Westcoast=[];
 RP250_Westcoast=[];
 RP500_Westcoast=[];
 
 
 for WestCoast=1:27;
 load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP1.mat']);
 RP1_Westcoast=[RP1_Westcoast;RP1]; 
 load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP2.mat']);
 RP2_Westcoast=[RP2_Westcoast;RP2];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP5.mat']);
 RP5_Westcoast=[RP5_Westcoast;RP5];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP10.mat']);
 RP10_Westcoast=[RP10_Westcoast;RP10];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP20.mat']);
 RP20_Westcoast=[RP20_Westcoast;RP20];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP25.mat']);
 RP25_Westcoast=[RP25_Westcoast;RP25];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP50.mat']);
 RP50_Westcoast=[RP50_Westcoast;RP50];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP100.mat']);
 RP100_Westcoast=[RP100_Westcoast;RP100];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP250.mat']);
 RP250_Westcoast=[RP250_Westcoast;RP250];
  load([DirOut filesep NAME_GRD_transects{WestCoast} '_RP500.mat']);
 RP500_Westcoast=[RP500_Westcoast;RP500];
 
 clear RP1 RP2 RP5 RP10 RP20 RP25 RP50 RP100 RP250 RP500
  
 
 end
 
  Name=[DirOut filesep 'RP1_Westcoast.mat'];
  save(Name,'RP1_Westcoast');
  Name=[DirOut filesep 'RP2_Westcoast.mat'];
  save(Name,'RP2_Westcoast');
  Name=[DirOut filesep 'RP5_Westcoast.mat'];
  save(Name,'RP5_Westcoast');
  Name=[DirOut filesep 'RP10_Westcoast.mat'];
  save(Name,'RP10_Westcoast');
  Name=[DirOut filesep 'RP20_Westcoast.mat'];
  save(Name,'RP20_Westcoast');
  Name=[DirOut filesep 'RP25_Westcoast.mat'];
  save(Name,'RP25_Westcoast');
  Name=[DirOut filesep 'RP50_Westcoast.mat'];
  save(Name,'RP50_Westcoast');
  Name=[DirOut filesep 'RP100_Westcoast.mat'];
  save(Name,'RP100_Westcoast');
  Name=[DirOut filesep 'RP250_Westcoast.mat'];
  save(Name,'RP250_Westcoast');
  Name=[DirOut filesep 'RP500_Westcoast.mat'];
  save(Name,'RP500_Westcoast');
  
  
  %Now can as a Tab delimited txt file
 headers={'StartLat','StartLon','EndLat','EndLon','TWL','TWLLat','TWLLon','DWL','DWLLat','DWLLon','zt','ztLat','ztLon','zc','zcLat','zcLon','MHW','MHWLat','MHWLon','Swash_Prob','Collision_Prob','OverWash_Prob','Inundation_Prob','DPY_Swash','DPY_Collision','DPY_OverWash','DPY_Inundation'};


%RP1
 for rr=1:length(headers);
     eval([headers{rr} '=RP1_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP1_WestCoast.txt']);
 
 %RP2
 for rr=1:length(headers);
     eval([headers{rr} '=RP2_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP2_WestCoast.txt']);
 
 %RP5
 for rr=1:length(headers);
     eval([headers{rr} '=RP5_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP5_WestCoast.txt']);
 
 %RP10
 for rr=1:length(headers);
     eval([headers{rr} '=RP10_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP10_WestCoast.txt']);
 
 %RP20
 for rr=1:length(headers);
     eval([headers{rr} '=RP20_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP20_WestCoast.txt']); 
 
 %RP25
 for rr=1:length(headers);
     eval([headers{rr} '=RP25_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP25_WestCoast.txt']);
 
 %RP50
 for rr=1:length(headers);
     eval([headers{rr} '=RP50_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP50_WestCoast.txt']); 
 
 %RP100
 for rr=1:length(headers);
     eval([headers{rr} '=RP100_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP100_WestCoast.txt']); 
 
 %RP250
 for rr=1:length(headers);
     eval([headers{rr} '=RP250_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP250_WestCoast.txt']);  
 
  %RP500
 for rr=1:length(headers);
     eval([headers{rr} '=RP500_Westcoast(:,' num2str(rr) ');']);
 end
 T=table(StartLat,StartLon,EndLat,EndLon,TWL,TWLLat,TWLLon,DWL,DWLLat,DWLLon,zt,ztLat,ztLon,zc,zcLat,zcLon,MHW,MHWLat,MHWLon,Swash_Prob,Collision_Prob,OverWash_Prob,Inundation_Prob,DPY_Swash,DPY_Collision,DPY_OverWash,DPY_Inundation);
 writetable(T,[DirOut filesep 'RP500_WestCoast.txt']);  
 
 