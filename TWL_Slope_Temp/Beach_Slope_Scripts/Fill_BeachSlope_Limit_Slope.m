%Code to fill in the NaN values for the Beach Slopes by looking up and
%down 0.5 km and taking an average of beach slopes of othee profiles 

% For use in updating the files within the StockSlope_Regional_v2 folder

%Load Slopes
Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional_v2'];
load([Dir filesep 'AllStockdon.mat']);

Prelim=AllStockdon.Slope;
CoName=AllStockdon.Counties;
Lat=AllStockdon.Y;
Lon=AllStockdon.X;
Prime=[];

Prelim(Prelim>0.1763)=NaN;
%Edit the Prelim to remove all slopes greater than threshold


NAME_GRD_transects  = {'San_Diego','Orange','Los_Angeles','Ventura','Santa_Barbara',...
    'Santa_Barbara_West','San_Luis_Obispo','Monterey','Santa_Cruz','San_Mateo','San_Francisco',...
    'Marin','Sonoma','Mendocino','Humboldt','Del_Norte','Curry','Coos','Douglas','Lane',...
    'Lincoln','Tillamook','Clatsop','Pacific','Grays_Harbor','Jefferson','Clallam'}; 


for ii = 1:length(Prelim);
    if ~isnan(Prelim(ii));
        Prime(ii)=Prelim(ii);
    else
        
        %so if it is indeed nan
        %First find the indices
        indx=[(ii-50):1:(ii+50)];
        
        indx(indx<0)=[];
        indx(indx>length(Prelim))=[];
        
        
        %Edit such that the neighboring profiles x and ymust be
        %within 500m of the middle point 
        indx(indx<1)=[];
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
        
        
        %If numbs is empty, then extend the window
        if isempty(numbs);
            indx=[(ii-100):1:(ii+100)];
        
        indx(indx<0)=[];
        indx(indx>length(Prelim))=[];

        tmpx=Lon(indx);
        tmpy=Lat(indx);
        
        orgx=Lon(ii);
        orgy=Lat(ii);
        
        %calculate the distance between and remove the appropriate indicies
        xcomp=(tmpx-orgx).^2;
        ycomp=(tmpy-orgy).^2;
        d=sqrt(xcomp+ycomp);
        
        a=find(d>1000);
        indx(a)=[];

        
        
        numbs=Prelim(indx);
        numbs(isnan(numbs))=[];
        numbs=mean(numbs);
        else
            numbs=mean(numbs);
        end

        
        Prime(ii)=numbs;
    end
end

        
AllStockdon.Slope_Limited=Prime';
save([Dir filesep 'AllStockdon.mat'],'AllStockdon');

