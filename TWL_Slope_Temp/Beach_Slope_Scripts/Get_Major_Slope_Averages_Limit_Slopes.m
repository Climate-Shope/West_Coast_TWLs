%Code to average the beach slope values at each major profile by looking up and
%down 50m and taking an average  

Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional_v2'];
load([Dir filesep 'AllStockdon.mat']);

Prelim=AllStockdon.Slope_Limited;
CoName=AllStockdon.Counties;
Lat=AllStockdon.Y;
Lon=AllStockdon.X;
Prime=[];

NAME_GRD_transects  = {'San_Diego','Orange','Los_Angeles','Ventura','Santa_Barbara',...
    'Santa_Barbara_West','San_Luis_Obispo','Monterey','Santa_Cruz','San_Mateo','San_Francisco',...
    'Marin','Sonoma','Mendocino','Humboldt','Del_Norte','Curry','Coos','Douglas','Lane',...
    'Lincoln','Tillamook','Clatsop','Pacific','Grays_Harbor','Jefferson','Clallam'}; 

for ii = 1:length(Prelim);
    if length(AllStockdon.Names{ii})>12;
        Prime(ii)=Prelim(ii);
    else
        
        %so if it is indeed nan
        %First find the indices
        indx=[(ii-5):1:(ii+5)];
        
        indx(indx<=0)=[];
        indx(indx>length(Prelim))=[];
        
        
        %Edit such that the neighboring profiles x and y must be
        %within 500m of the middle point 
        tmpx=Lon(indx);
        tmpy=Lat(indx);
        
        orgx=Lon(ii);
        orgy=Lat(ii);
        
        %calculate the distance between and remove the appropriate indicies
        xcomp=(tmpx-orgx).^2;
        ycomp=(tmpy-orgy).^2;
        d=sqrt(xcomp+ycomp);
        
        a=find(d>500);
        indx(a)=[];
        
        numbs=Prelim(indx);
        numbs(isnan(numbs))=[];
        numbs=mean(numbs);
        Prime(ii)=numbs;
        
        
        

        Prime(ii)=numbs;
    end
end

        
AllStockdon.Slope_MeanMajors_Limited=Prime';
save([Dir filesep 'AllStockdon.mat'],'AllStockdon');

