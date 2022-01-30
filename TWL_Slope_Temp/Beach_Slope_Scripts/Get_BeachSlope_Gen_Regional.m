%Load the wave conditions for a specific profile to extract the maximum
%energy conditions.

% This is a second variation of the calculated beach slope and shoudl be
% run after v5

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

%Setup for the process of getting the slopes for intput into the Stockdon et al. (2006) runup method for each
%location. This recoreds the slopes as two separate
%vectors, one for the names of each minor transect and the actual transect,
%then another for their values. 

% Create the out directory 
DirOut=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional_v2'];
if ~exist(DirOut,'dir'), mkdir(DirOut); end

%Do not need NTRs for this section 

%Initialize vectors
AllSlopes=[];
AllNames={};
AllCounties={};
AllLat=[];
AllLon=[];


for WestCoast=1:length(NAME_GRD_transects);

%% Determine how many files to loop through
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\time_series_reconstructed',filesep,'ii_grid_array.mat']);
Ntr=ii_grid_array;

load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\Profile_ESI', filesep, NAME_GRD_transects{WestCoast}, '_ShoreType.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\Profile_ESI',filesep,NAME_GRD_transects{WestCoast}, '_SeaIntersect.mat']);
load(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects{WestCoast}, '\Profile_ESI',filesep,NAME_GRD_transects{WestCoast}, '_ESI_Array.mat']);

utmzone=UTM{WestCoast};

%load in the filler slopes determined from regional (such as southern CA)
%areas
load('F:\West_Coast_TWL_Hazards\01_Data\CA_Slope_Data\Fill_Slopes.mat');
cd('F:\West_Coast_TWL_Hazards\_STEP\TWL_Slope_Temp');
SlopeArray=[];

for ii=ii_grid_array;
    OBJECTID=dec2base(ii,10,4);
    if ii==1215;
        pause=1;%Just for testing specific profiles
    end
   

%If a profile exists where the toe locaiton had special fixes load that profile, if not, load the normal profile    
if exist(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'Morphology_Fix_Toe', filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat'],'file');
    fname=['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'Morphology_Fix_Toe', filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat'];
else
    fname=['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'Morphology', filesep, 'SlopeInfo_Transect_', num2str(ii), '.mat'];
end
if exist(fname,'file');
    load(fname);    
    prof_name=['profile_', dec2base(ii,10,4)];
    
    
%%    %detect whether the profile is in deg or UTM
    if ~isempty(SlopeInfo.ypf);
        tam=round(SlopeInfo.ypf(1));
        tam=numel(num2str(tam));
        if tam<=3;
            [x,y,utmz] = deg2utm(SlopeInfo.ypf,SlopeInfo.xpf);
            SlopeInfo.ypf=y;
            SlopeInfo.xpf=x;
        end
    end
    
%%    %Grab the shoreline type (ESI)
    ShoreID=SlopeInfo.ESI;
    if isempty(ShoreID);
        ESIID={NaN};
    else
        ESIID=ShoreID;
    end
    
    ESI=ESIID;

% Flag for seawall
   Rey=SlopeInfo.Rey;

%% Slope Analysis     
    if ~isempty(SlopeInfo.depth) && SlopeInfo.maxi~=0;
     %Load the first point in the profile 
    x=SlopeInfo.xpf(1);
    y=SlopeInfo.ypf(1);
    
    [Lat,Lon] = utm2deg(x,y,utmzone);
    
    %Pick a fill_slope from reigonal values
    xcomp=(Lat-Fill_Slopes(:,1)).^2;
    ycomp=(Lon-Fill_Slopes(:,2)).^2;
    d=sqrt(xcomp+ycomp);
    [~,I]=min(abs(d));    
    fslope=Fill_Slopes(I,:);
    
    
    %Determine Slope and Runup formulation
     if max(SlopeInfo.depth)<3
        StockTemp=NaN;
     else

        StockTemp=SlopeInfo.Runup_Method.StockSlope_reg;
     end
    else
        StockTemp=NaN;
    end
    
else
        
      StockTemp=NaN;  
      ShoreID=find(ShoreType.OBJECTID==ii);
    if isempty(ShoreID);
        ESIID={NaN};
    else
        if isnan(ShoreType.ESI{ShoreID,1});
            ESIID={NaN};
        else
        ESIID=cellstr(ShoreType.ESI{ShoreID,1});
        end
    end
       prof_name=['profile_', dec2base(ii,10,4)];
       Rey=[];

end
    
    
   %Get X and Y points 
   
   if isempty(SlopeInfo.xpf);
   tmpx=[NaN];
   tmpy=[NaN];
   else
   tmpx=[SlopeInfo.xpf(1)];
   tmpy=[SlopeInfo.ypf(1)];
   end 
    
%% Now loop through the subprofiles and get their slopes as well
    %If the profle is 1, then do the first 5 in the folder
    %If at the end, do the last 5
    %Otherwise, look for the files that are 5 before and 5 after
    minordir=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\interp_minor_profiles_wave_model_edited'];
    zip=[5:9,1:5];
    Rey=[];
    if ii==1;
        vsl=6:10;
    elseif ii==ii_grid_array(end);
        vsl=2:5;
    else
        vsl=2:10;
    end
    minor_prof_name={};
    for jj=vsl;
        if jj<6;
            minorname=[minordir,filesep,'profile_', dec2base(ii-1,10,4),'_',num2str(zip(jj)),'.mat'];
            clap=['profile_', dec2base(ii-1,10,4),'_',num2str(zip(jj))];
        else
            minorname=[minordir,filesep,'profile_', dec2base(ii,10,4),'_',num2str(zip(jj)),'.mat'];
            clap=['profile_', dec2base(ii,10,4),'_',num2str(zip(jj))];
        end
        
        if ~exist(minorname,'file');
            if jj<6;
                oldname=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\interp_minor_profiles_wave_model',filesep,'profile_', dec2base(ii-1,10,4),'_',num2str(zip(jj)),'.mat'];
                if exist(oldname,'file');
                copyfile(oldname,minordir)
                end
          
            else
                oldname=['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects{WestCoast} '\interp_minor_profiles_wave_model',filesep,'profile_', dec2base(ii,10,4),'_',num2str(zip(jj)),'.mat'];
                if exist(oldname,'file');
                copyfile(oldname,minordir)
                end
            end
        end
        
        if exist(minorname,'file');
            load(minorname);
        
            
            %If Santa Barbara West region, need to convert xy to UTM
            if strcmp(NAME_GRD_transects{WestCoast},'Santa_Barbara_West');
                [wave_profile.xpf,wave_profile.ypf,wave_profile.utmzone] = deg2utm(wave_profile.ypf,wave_profile.xpf);
            end
        
        if jj<6;
            if isempty(wave_profile.depth) || length(wave_profile.depth)<30
                copyfile(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'interp_minor_profiles_wave_model', filesep, 'profile_', dec2base(ii-1,10,4),'_',num2str(zip(jj)),'.mat'],minordir)
                load(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'interp_minor_profiles_wave_model', filesep, 'profile_', dec2base(ii-1,10,4),'_',num2str(zip(jj)),'.mat']);
            end
        else
            if isempty(wave_profile.depth) || length(wave_profile.depth)<30
                copyfile(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'interp_minor_profiles_wave_model', filesep, 'profile_', dec2base(ii,10,4),'_',num2str(zip(jj)),'.mat'],minordir)
                load(['F:\West_Coast_TWL_Hazards\03_Results' filesep NAME_GRD_transects{WestCoast} filesep, 'interp_minor_profiles_wave_model', filesep, 'profile_', dec2base(ii,10,4),'_',num2str(zip(jj)),'.mat']);
            end
        end
        
        minor_prof_name=[minor_prof_name;clap];
        
        
    if ~isempty(wave_profile.ypf);
        tam=round(wave_profile.ypf(1));
        tam=numel(num2str(tam));
        if tam<=3;
            [x,y,utmz] = deg2utm(wave_profile.ypf,wave_profile.xpf);
            wave_profile.ypf=y';
            wave_profile.xpf=x';
        end
    end
%% Now clip the profile to 300m onshore 
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
        
        if isempty(wave_profile.depth);
            StockTemp=[StockTemp;NaN];
            tmpx=[tmpx;NaN];
            tmpy=[tmpy;NaN];
            continue
        end
        
        if max(wave_profile.depth)<3
            StockTemp=[StockTemp;NaN];
        else
            
            
%% Calculate the critical thresholds of the elevation profile 
         [SlopeInfo]=MorphLove_RegSlope(wave_profile,ESI,ESI_Array,Rey,OBJECTID,UTM{WestCoast});
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
        try
           [SlopeInfo]=MorphLove_RegSlope(wave_profile,ESI,ESI_Array,Rey,OBJECTID,UTM{WestCoast});
        catch
            flag=1;%If the edited profile fails
        end
    end

%% If the flag for an erroneous profile does not exist, calculate the slope
%  otherwise, NaN
    if ~exist('flag','var')
        wave_profile.depth=SlopeInfo.depth;
        wave_profile.xpf=SlopeInfo.xpf;
        wave_profile.ypf=SlopeInfo.ypf;
        wave_profile.L=SlopeInfo.L;
    [RunupMethod]=Runup_Method(SlopeInfo,wave_profile,Rey,fslope);
        StockTemp=[StockTemp;RunupMethod.StockSlope_reg];
    else
         StockTemp=[StockTemp;NaN];
         clear flag
    end
        
        
        end
        else
         minor_prof_name=[minor_prof_name;clap];
         StockTemp=[StockTemp;NaN];
        end
        
        if isempty(wave_profile.xpf);
        tmpx=[tmpx;NaN];
        tmpy=[tmpy;NaN];
        else
        tmpx=[tmpx;wave_profile.xpf(1)];
        tmpy=[tmpy;wave_profile.ypf(1)];
        end
        
    end
    StockTemp(StockTemp==fslope(:,3))=NaN;
    
    names=[prof_name;minor_prof_name];
    
    if length(StockTemp)==1 && isnan(StockTemp);
        StockTemp=nan(length(names),1);
        tmp3=repmat(NAME_GRD_transects{WestCoast},length(names),1);
        tmp3=cellstr(tmp3);
    else
    %Reorder based on the ii and vsl
    if ii==1 %order doesnt need to change
        %names={prof_name;
    elseif ii==ii_grid_array(end)%first value needs to be placed at the end
        StockTemp=[StockTemp;StockTemp(1)];
        StockTemp(1)=[];
        names=[names;names{1}];
        names(1)=[];
        tmpx=[tmpx;tmpx(1)];
        tmpx(1)=[];
        tmpy=[tmpy;tmpy(1)];
        tmpy(1)=[];
        
    else %first value needs to be placed in the middle
        tmp1=StockTemp(2:5);
        tmp2=StockTemp(6:end);
        StockTemp=[tmp1;StockTemp(1);tmp2];
        
        tmp1=names(2:5);
        tmp2=names(6:end);
        names=[tmp1;names{1};tmp2];
        
        tmp1=tmpx(2:5);
        tmp2=tmpx(6:end);
        tmpx=[tmp1;tmpx(1);tmp2];
        
        tmp1=tmpy(2:5);
        tmp2=tmpy(6:end);
        tmpy=[tmp1;tmpy(1);tmp2];
              
    end
    
    %now need the counties
    tmp3=repmat(NAME_GRD_transects{WestCoast},length(StockTemp),1);
    tmp3=cellstr(tmp3);

    end
    
%% Combine into continuous vectors across all counties       
AllSlopes=[AllSlopes;StockTemp];
AllNames=[AllNames;names];
AllCounties=[AllCounties;tmp3];   
AllLat=[AllLat;tmpy];
AllLon=[AllLon;tmpx];    
    
    
%% Catch for some error values       
AllSlopes(AllSlopes<0)=NaN;
AllSlopes(AllSlopes>tan(6.3*(pi/180)))=NaN;

end


%% Create Output variable 
AllStockdon.Slope=AllSlopes;
AllStockdon.Names=AllNames;
AllStockdon.Counties=AllCounties;
AllStockdon.X=AllLon;
AllStockdon.Y=AllLat;

save([DirOut filesep NAME_GRD_transects{WestCoast} '_AllStockdon.mat'],'AllStockdon'); 

AllSlopes=[];
AllNames={};
AllCounties={};
AllLat=[];
AllLon=[];
clear AllStockdon

end




    