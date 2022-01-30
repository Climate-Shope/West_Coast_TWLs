%Code to average the Stockdon beach slope values at each major profile by looking up and
%down 50m and taking an average of whatever turns up. 

Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional'];
load([Dir filesep 'AllStockdon.mat']);

Prelim=AllStockdon.Slope_Edited;
CoName=AllStockdon.Counties;
Lat=AllStockdon.Y;
Lon=AllStockdon.X;
Prime=[];

NAME_GRD_transects  = {'San_Diego','Orange','Los_Angeles','Ventura','Santa_Barbara',...
    'Santa_Barbara_West','San_Luis_Obispo','Monterey','Santa_Cruz','San_Mateo','San_Francisco',...
    'Marin','Sonoma','Humboldt','Del_Norte''Del_Norte','Curry','Coos','Douglas','Lane',...
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
        
        
        %OKay so now edit such that the neighboring profiles x and ymust be
        %within 500m of the middle point Need to build in a fix for SB
        %transition, may skip Santa Barbara West and leave it with the
        %normal indicies.
        if strcmp(CoName{ii},'Santa_Barbara') || strcmp(CoName{ii},'Santa_Barbara_West')
        else
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
        end
        
        numbs=Prelim(indx);
        numbs(isnan(numbs))=[];
        numbs=mean(numbs);
        Prime(ii)=numbs;
        %numbs should not be empty by this new schematization 
        
 
        
        Prime(ii)=numbs;
    end
end

        
AllStockdon.Slope_MeanMajors=Prime';
save([Dir filesep 'AllStockdon.mat'],'AllStockdon');

