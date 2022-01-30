function [Out]=Runup_Method_Update(SlopeInfo,wave_profile,fslope);
%function to select the correct runup method along each profile 
%%%% List of potential points
%maxi = location of the maxima closest to shore
maxi=SlopeInfo.maxi;
maxi_ele=SlopeInfo.maxi_ele;

%toe_loc = location of toe of the strcutre cliff, dune based on what the
%water interacts with from stockdon 
toe_loc=SlopeInfo.toe_loc;
toe_ele=SlopeInfo.toe_ele;

if toe_loc<1;
    toe_loc=1;
    toe_ele=wave_profile.depth(1);
end


%slope_break_locs = locations of the slope breaks in the profile, need to
%edit them to figure out if they are relevant

CliffTop_Loc=SlopeInfo.CliffTop_Loc;
CliffTop_Ele=SlopeInfo.CliffTop_Ele;
CliffJun_Loc=SlopeInfo.CliffJun_Loc;
CliffJun_Ele=SlopeInfo.CliffJun_Ele;

berm_loc=SlopeInfo.berm_loc;%If there is a rocky berm/bench

%first need to figure out the distance between the MHW mark and toe (if 
%extant)

MHW=SlopeInfo.MHW;
MHW_loc=SlopeInfo.MHW_loc;


% Now ESI
ESI=SlopeInfo.ESI;

% enter a flag for detached seawall
Rey=SlopeInfo.Rey;
if ~isempty(Rey);
    SeaWall=1;
else
    SeaWall=[];
end

%% Stocktoe is the lower of the two potential toes and is used for an inital 
%  slope calculation for the Stockdon Runup calculation
if SlopeInfo.toe_onshore >= 1 && SlopeInfo.toe_loc >= 1

if SlopeInfo.depth(SlopeInfo.toe_onshore)< SlopeInfo.depth(SlopeInfo.toe_loc)
    Stocktoe=SlopeInfo.toe_onshore;
elseif SlopeInfo.depth(SlopeInfo.toe_onshore)> SlopeInfo.depth(SlopeInfo.toe_loc)
    Stocktoe=SlopeInfo.toe_loc;
else
    Stocktoe=SlopeInfo.toe_onshore;
end
else
    Stocktoe=-1;
end


%% First, Calculate the local stockdon slope if applicable

%in simplest scenario it will just be MHW to the Stocktoe variable

if Stocktoe > MHW_loc && Stocktoe > 1 && Stocktoe-MHW > 5;
    StockSlope=(SlopeInfo.depth(Stocktoe)-MHW)/(Stocktoe-MHW_loc);
elseif Stocktoe > MHW_loc && Stocktoe > 1 && Stocktoe-MHW < 5;
    StockSlope=(SlopeInfo.depth(Stocktoe)-SlopeInfo.depth(1))/(Stocktoe-1);
elseif Stocktoe < 1; %Here want to go from MHW to the maxi or cliffjun
    if ~isempty(SlopeInfo.maxi) && ~isempty(SlopeInfo.CliffJun_Loc);
        if SlopeInfo.maxi<SlopeInfo.CliffJun_Loc
            loc=SlopeInfo.maxi;
        else
            loc=SlopeInfo.CliffJun_Loc;
        end
        StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    elseif isempty(SlopeInfo.maxi) && ~isempty(SlopeInfo.CliffJun_Loc);
            loc=SlopeInfo.CliffJun_Loc;
            StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    elseif ~isempty(SlopeInfo.maxi) && isempty(SlopeInfo.CliffJun_Loc);
            loc=SlopeInfo.maxi;
             StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    end
elseif Stocktoe==MHW_loc;%If the toe and the MHW location ar coincident, usually menaing no actual toe could be found
    if ~isempty(SlopeInfo.maxi) && ~isempty(SlopeInfo.CliffJun_Loc);
        if SlopeInfo.maxi<SlopeInfo.CliffJun_Loc
            loc=SlopeInfo.maxi;
        else
            loc=SlopeInfo.CliffJun_Loc;
        end
        StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    elseif isempty(SlopeInfo.maxi) && ~isempty(SlopeInfo.CliffJun_Loc);
            loc=SlopeInfo.CliffJun_Loc;
            StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    elseif ~isempty(SlopeInfo.maxi) && isempty(SlopeInfo.CliffJun_Loc);
            loc=SlopeInfo.maxi;
             StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    end
elseif Stocktoe<MHW_loc; %in this case, use the slope from MHW to the junciton or maximum because the toe should not be seaward of MHW
    if ~isempty(SlopeInfo.maxi) && ~isempty(SlopeInfo.CliffJun_Loc);
        if SlopeInfo.maxi<SlopeInfo.CliffJun_Loc
            loc=SlopeInfo.maxi;
        else
            loc=SlopeInfo.CliffJun_Loc;
        end
        StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    elseif isempty(SlopeInfo.maxi) && ~isempty(SlopeInfo.CliffJun_Loc);
            loc=SlopeInfo.CliffJun_Loc;
            StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    elseif ~isempty(SlopeInfo.maxi) && isempty(SlopeInfo.CliffJun_Loc);
            loc=SlopeInfo.maxi;
             StockSlope=(SlopeInfo.depth(loc)-MHW)/(loc-MHW_loc);
    end
end


%Catch to determine if the slope should be used in the reigonal
%calculation 

if (atan(StockSlope)*(180/pi))<34%threshold determined by testing and fall around the angle of repose for medium sand
    StockSlope_reg=StockSlope;
else
    StockSlope_reg=NaN;
end

%% Get a few parameters out of the way
if ~exist('StockSlope','var');
    StockSlope=NaN;
end




%% Now to add in reduction factors and methodologies

%Roughness 1st, these gamma values for TAW are from Allan and others 2015
%from their oregon report and the FEMA guidelines for coastal flood hazard
%analysis (NHC 2005).
if strcmp(ESI,'1A') || strcmp(ESI,'1B')|| strcmp(ESI,'2A')|| strcmp(ESI,'8A')|| strcmp(ESI,'8B')||strcmp(ESI,'8F')
    gamma=1.0;
    Method='TAW';
elseif strcmp(ESI,'6D') || strcmp(ESI,'6B') || strcmp(ESI,'6A')
    if strcmp(ESI,'6D')
      gamma=0.65;  
    elseif strcmp(ESI,'6B') || strcmp(ESI,'8C')
      gamma=0.55;
    else
        gamma=0.7;
    end
    Method='TAW';
else
    Method='Stockdon';
    gamma=NaN;
end



%% fix up the initial slope
if (atan(StockSlope)*(180/pi))>34;%A coarse cutoff for where TAW should be used versus Stockdon based on angle of repose
    Method='TAW'; 
else
    Method='Stockdon';
end


%% Define if the backshore has a different ESI type (such as a beach backed by a cliff)



if exist('CliffJun_Ele','var') && CliffJun_Ele~=0;
   temp=CliffJun_Ele-SlopeInfo.toe_ele;
   if temp > 10;
      Backshore_ESI='1A';
   else
      Backshore_ESI='0';
   end
else
     Backshore_ESI='0';
end
%define gamma for the backing cliff/bluff. Assumed gamma of just a smooth
%cliff
if exist('Backshore_ESI','var')
    if strcmp(Backshore_ESI,'1A')
    backgamma=1.0;
    ExceedMethod='TAW';
    else
    backgamma=NaN;
    end
else
    backgamma=NaN;
end

%% Next up, calc the slope around the further inland toe to determine the exceed method
%Do this via 2 methods
%1 slope of the toe and the next 20 points to see if there is a
%steep part (determined by testing)
if toe_loc+20 < length(wave_profile.depth);
    slope_range=wave_profile.depth(toe_loc:toe_loc+20);
else
    slope_range=wave_profile.depth(toe_loc:end);
end
slope_range=diff(slope_range);
a=find(slope_range>0.67);
if ~isempty(a);
    one='TAW';
else
    one='Stockdon';
end


%now to do the slope calc between the toe and the cliffjunction 
if wave_profile.depth(SlopeInfo.CliffJun_Loc)-wave_profile.depth(toe_loc) > 10
    slope=(wave_profile.depth(SlopeInfo.CliffJun_Loc)-wave_profile.depth(toe_loc))/(SlopeInfo.CliffJun_Loc-toe_loc);
    if slope>0.67
       two='TAW';
    else
        two='Stockdon';
    end 
else
    two='Stockdon';
end
%If both methods read TAW, the runup method once the toe is exceeded is TAW
if strcmp(one,'TAW') || strcmp(two,'TAW')
    ExceedMethod='TAW';
else
    ExceedMethod='Stockdon';
end

%Need to modify this a bit for specific things like 6B and 5
if strcmp(ESI,'6B') || ~isempty(Rey)
    ExceedMethod='TAW';
end

if strcmp(ESI,'2A') || ~isempty(Rey)
    ExceedMethod='TAW';
end

if strcmp(ESI,'1B') || ~isempty(Rey)
    ExceedMethod='TAW';
end


%% Now want to check for instaces where there is a stockdon slope defined, but poorly defined because of the short length
%  and redefine the stockdon slope based on initial input
if strcmp(Method,'Stockdon') && strcmp(ESI,'6B')
    if SlopeInfo.toe_loc<10
        StockSlope=fslope(:,3);
    end
end

if SlopeInfo.toe_loc<1
    Method='TAW';
end

%If normal method is TAW, then exceed method must be TAW by definition
if strcmp(Method,'TAW');
    ExceedMethod='TAW';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Next add a piece to help deal with Seawalls that are so low lying that it probably needs to 
% be just stockdon 

if ~isempty(SeaWall) && strcmp(ExceedMethod,'TAW');
    %Grab the slope near the toe and the overtopping points
    if SlopeInfo.toe_loc<1
    tt=(SlopeInfo.depth(SlopeInfo.overtop_point)-SlopeInfo.depth(1))/(SlopeInfo.overtop_point-1);
    tt=atan(tt)*(180/pi);        
    else
    tt=(SlopeInfo.depth(SlopeInfo.overtop_point)-SlopeInfo.depth(SlopeInfo.toe_loc))/(SlopeInfo.overtop_point-SlopeInfo.toe_loc);
    tt=atan(tt)*(180/pi);
    end
    if tt<10
       ExceedMethod='Stockdon';  
    end
    
end

%% Produce output variable
Out.Method=Method; 
if exist('ExceedMethod','var')
    Out.ExceedMethod=ExceedMethod;
else
    Out.ExceedMethod=Method;
end
Out.StockSlope=StockSlope;
if SlopeInfo.toe_loc<1 
%If the storckdon toe is less than one (menaing it does not exist, use a
%regional estimate
   Out.StockSlope=fslope(:,3);
end

%%% Now add a parameter to edit the stockslope such that if the toe used to calculate the slope is less than 10, it sets it to the reigonal slope 
if Stocktoe<10;
     Out.StockSlope=fslope(:,3);
end

Out.StockSlope_reg=StockSlope_reg;
Out.gamma=gamma;
Out.backgamma=backgamma;
Out.ESI=ESI;
Out.BackshoreESI=Backshore_ESI;
Out.SeaWall=SeaWall;

end






















