%Generate an updated form of the tranects from hand edited shapefile
%generated by prior steps and corrected in ArcGIS
clear all; run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 
NAME_GRD_transects  = 'Douglas'; 
run('run_init_directories.m')

transect_dir='E:\West_Coast_TWL_Hazards\01_Data\Transects';

S=shaperead(['E:\West_Coast_TWL_Hazards\_STEP\_' NAME_GRD_transects filesep NAME_GRD_transects '_Transects' filesep  NAME_GRD_transects '_Transects_v3.shp']);


X=[];
Y=[];

%Get start and end transect points
for ii=1:numel(S);
    x=S(ii).X;
    y=S(ii).Y;

    xvec=[x(1) x(end-1)];
    yvec=[y(1) y(end-1)];
    
    X=[X;xvec];
    Y=[Y;yvec];

end


%add length onshore
xx=[];
yy=[];
for ii=1:length(X(:,1));
    u=X(ii,2)-X(ii,1);
    v=Y(ii,2)-Y(ii,1);
    
    [r,az] = pcoord(u,v);
    %Adding 500 m onshore from the shoreline point
    r=500;
    
    [u,v] = rcoord(r, az);
    x=X(ii,2)+u;
    y=Y(ii,2)+v;

    xx=[xx;x];
    yy=[yy;y];


end



%add length offshore
for ii=1:length(X(:,1));
    u=X(ii,2)-X(ii,1);
    v=Y(ii,2)-Y(ii,1);
    
    [r,az] = pcoord(u,v);
    %Add 1500 m offshore from shoreline point
    r=1500;
    
    [u,v] = rcoord(r, az);
    x=X(ii,1)-u;
    y=Y(ii,1)-v;

    X(ii)=x;
     Y(ii)=y;


end




X=[X xx];
Y=[Y yy];

name=[NAME_GRD_transects '_Transects_v3'];
eval([name '.X=X;']);
eval([name '.Y=Y;']);

save([transect_dir filesep name],name);
    









