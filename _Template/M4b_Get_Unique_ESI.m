%Edits and saves on value of ESI and armoring data for a transect if multiple options are
%possible, primarily, which one that the transect intersects first from
%offshore to oonshore

warning off all 

NAME_GRD_transects  = 'Douglas'; 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');

file_esi          = ['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects, '_ESI_Array.mat'];
load(file_esi);

load(['F:\West_Coast_TWL_Hazards\_STEP\_', NAME_GRD_transects,filesep, 'Transects2.mat']);


temp=size(ESI_Array);
temp=temp(1);

OBJECTID=[];
ESIIS=[];
X=[];
Y=[];
for ii=1:temp;
    
    a=ESI_Array{ii,1};
    OBJECTID=[OBJECTID;a];
    a=ESI_Array{ii,2};
    X=[X;a];
    a=ESI_Array{ii,3};
    Y=[Y;a];    
    a=ESI_Array{ii,4};
    if ~isnan(a)
    a=cellstr(a);
    ESIIS=[ESIIS;a];  
    else
    ESIIS=[ESIIS;{NaN}];  
    end
end

temp=unique(OBJECTID);

ESIIS_2=[];
for ii=1:length(temp);
    
a=find(OBJECTID==ii);
if isempty(a)
    ESIIS_2=[ESIIS_2;{NaN}];
else
if length(a)==1;
    b=ESIIS{a,1};
    if ~isnan(b)
    b=cellstr(b);
    ESIIS_2=[ESIIS_2;b];
    else
    ESIIS_2=[ESIIS_2;{NaN}];
    end
    
else
    %Want the point closest to the transect head
    tX=X(a);
    tY=Y(a);
    xoff=Transects2.xpf(ii,1);
    yoff=Transects2.ypf(ii,1);
    
    xcomp=(tX-xoff).^2;
    ycomp=(tY-yoff).^2;
    distance=sqrt(xcomp+ycomp);
    [~,I]=min(distance);
    a=a(I);
    b=ESIIS{a,1};
    b=cellstr(b);
    ESIIS_2=[ESIIS_2;b];
end


end
end

ShoreType.ESI=ESIIS_2;
ShoreType.OBJECTID=temp;

save(['F:\West_Coast_TWL_Hazards\03_Results\', NAME_GRD_transects, '\Profile_ESI',filesep,NAME_GRD_transects,'_ShoreType'],'ShoreType');  











%        
       
       