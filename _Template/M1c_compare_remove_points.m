%Converts Transects from M1b to generate a non-repeated and complete matlab
%file verison of transect locations
clear all; run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Douglas'; 
run('run_init_directories.m')

transect_dir='E:\West_Coast_TWL_Hazards\01_Data\Transects';
point_dir='E:\West_Coast_TWL_Hazards\01_Data\Transects\Transect_Center_Points';

%load tranects
load([transect_dir filesep NAME_GRD_transects '_Transects.mat']);
eval(['transects=' NAME_GRD_transects '_Transects;']);
%load points
S=shaperead([point_dir filesep NAME_GRD_transects '_Centers.shp']);

X=[];
Y=[];
for ii=1:numel(S);
    x=round(S(ii).F1);
    y=round(S(ii).F2);
    
    X=[X;x];
    Y=[Y;y];
end

Point_Rows=[X,Y];
Trans_Rows=[round(transects.X(:,2)), round(transects.Y(:,2))];

idx=ismember(Trans_Rows,Point_Rows,'rows','legacy');

z=find(idx);

VIS=0;
if VIS==1;
figure;hold on
scatter(Trans_Rows(:,1),Trans_Rows(:,2));
scatter(Point_Rows(:,1),Point_Rows(:,2));

end

transects2=transects;

transects2.X=transects2.X(z,:);
transects2.Y=transects2.Y(z,:);

name=[NAME_GRD_transects '_Transects_v2'];

eval([name '=transects2;']);

save([transect_dir filesep name,name]);
