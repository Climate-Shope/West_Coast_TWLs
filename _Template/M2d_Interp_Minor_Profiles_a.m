%Same as M2C but calculates profilie elevations for 10-m spaced minor
%profiles

%Part 1 casts the transects
%Part 2 does the interpolations 
warning off all 
NAME_GRD_transects  = 'Douglas'; 
addpath('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
dirout      = ['E:\West_Coast_TWL_Hazards\_STEP\_' NAME_GRD_transects];


file_MHW                =['E:\West_Coast_TWL_Hazards\03_Results\MHW' filesep NAME_GRD_transects '_MHW.mat']; 
file_transects          = ['E:\West_Coast_TWL_Hazards\_STEP\_' NAME_GRD_transects filesep 'Transects2.mat'];


num_dem=1;



h1 = 0; 
h2 = +55; % profile limits, bathymetry <0 
resolution_m = 1; % m 
prof_cut_len=50;
%% Step 1, loadin in MHW and calc distance between shorleine points. 

load(file_MHW);

%calculate euclidean distance and assume any points that are greater than
%150 m apart represent an intentional omission for harbors, river mouths,
%etc. 


xcomp=diff(MHW(:,1)).^2;
ycomp=diff(MHW(:,2)).^2;
distance=sqrt(xcomp+ycomp);

skip=find(distance>150);


%% Step 2, cast the transects

load(file_transects);
tempx=[];
tempy=[];
for ii=1:length(Transects2.xpf(:,1))-1;
    if ismember(ii,skip);
        continue
    end
    
    offx=linspace(Transects2.xpf(ii,1),Transects2.xpf(ii+1,1),11);
    offx([1;11])=[];
    offy=linspace(Transects2.ypf(ii,1),Transects2.ypf(ii+1,1),11);
    offy([1;11])=[];
    mnr=1:length(offy);


    onx=linspace(Transects2.xpf(ii,2),Transects2.xpf(ii+1,2),11);
    onx([1;11])=[];
    ony=linspace(Transects2.ypf(ii,2),Transects2.ypf(ii+1,2),11);
    ony([1;11])=[];


    maj=repmat(ii,length(onx),1);
    
    
    outx=[offx',onx',maj,mnr'];
    outy=[offy',ony',maj,mnr'];
    
    tempx=[tempx;outx];
    tempy=[tempy;outy];
end


%% Step 3, extend transects on and offshore

xx=[];
yy=[];
for ii=1:length(tempx(:,1));
    u=tempx(ii,2)-tempx(ii,1);
    v=tempy(ii,2)-tempy(ii,1);
    
    [r,az] = pcoord(u,v);
    r=100;
    
    [u,v] = rcoord(r, az);
    x=tempx(ii,2)+u;
    y=tempy(ii,2)+v;

    xx=[xx;x];
    yy=[yy;y];


end

tempx(:,6)=xx;
tempy(:,6)=yy;


%add length offshore
xx=[];
yy=[];
for ii=1:length(tempx(:,1));
    u=tempx(ii,2)-tempx(ii,1);
    v=tempy(ii,2)-tempy(ii,1);
    
    [r,az] = pcoord(u,v);
    r=100;
    
    [u,v] = rcoord(r, az);
    x=tempx(ii,1)-u;
    y=tempy(ii,1)-v;

    xx=[xx;x];
    yy=[yy;y];

end

tempx(:,5)=xx;
tempy(:,5)=yy;



%% Now Save as a new Transect3 Listing

Transects3.xpf=tempx;
Transects3.ypf=tempy;

save([dirout filesep 'Transects3'],'Transects3');



transects_mod.xoff = Transects3.xpf(:,5);
transects_mod.xon = Transects3.xpf(:,6);
transects_mod.yoff = Transects3.ypf(:,5);
transects_mod.yon = Transects3.ypf(:,6);
transects_mod.MAJORID=Transects3.xpf(:,3);
transects_mod.MINORID=Transects3.xpf(:,4);
INCORRECT_RUNS = [];


    Transects4.xpf=[transects_mod.xoff transects_mod.xon];
    Transects4.ypf=[transects_mod.yoff transects_mod.yon];
    Transects4.OBJECTID=transects_mod.MAJORID;
    Transects4.MINORID=transects_mod.MINORID;



save('Transects4','Transects4');

























