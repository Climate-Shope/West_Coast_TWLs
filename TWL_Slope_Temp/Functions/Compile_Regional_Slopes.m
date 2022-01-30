%code to append the slopes by county into one large vector

Dir=['F:\West_Coast_TWL_Hazards\03_Results\StockSlope_Regional_v2'];

NAME_GRD_transects  = {'San_Diego','Orange','Los_Angeles','Ventura','Santa_Barbara',...
    'Santa_Barbara_West','San_Luis_Obispo','Monterey','Santa_Cruz','San_Mateo','San_Francisco',...
    'Marin','Sonoma','Mendocino','Humboldt','Del_Norte','Del_Norte','Curry','Coos','Douglas','Lane',...
    'Lincoln','Tillamook','Clatsop','Pacific','Grays_Harbor','Jefferson','Clallam'};

Name=[];
CoName=[];
Slope=[];
Lat=[];
Lon=[];

for ii=1:length(NAME_GRD_transects);
    load([Dir filesep NAME_GRD_transects{ii} '_AllStockdon.mat']);
    
    if AllStockdon.Y(1)>AllStockdon.Y(end);
        AllStockdon.Names=flipud(AllStockdon.Names);
        AllStockdon.Counties=flipud(AllStockdon.Counties);
        AllStockdon.Slope=flipud(AllStockdon.Slope);
        AllStockdon.Y=flipud(AllStockdon.Y);
        AllStockdon.X=flipud(AllStockdon.X);
    end
        
    Name=[Name;AllStockdon.Names];
    CoName=[CoName;AllStockdon.Counties];
    Slope=[Slope;AllStockdon.Slope];
    Lat=[Lat;AllStockdon.Y];
    Lon=[Lon;AllStockdon.X];

end


AllStockdon.Slope=Slope;
AllStockdon.Names=Name;
AllStockdon.Counties=CoName;
AllStockdon.X=Lon;
AllStockdon.Y=Lat;

save([Dir filesep 'AllStockdon.mat'],'AllStockdon');




