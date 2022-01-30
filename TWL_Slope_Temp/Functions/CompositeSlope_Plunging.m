function [Out] = CompositeSlope_Plunging(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);

%script derived from code utilized by Allan and others (2015)

% code to calculate the composite slope in use for TAW, the function 
% utilizes the greatest wave condition at each transect to define this slope

%In this case, the slope will be set as the local slope defined by the
%stockdon TWL and the Stockdon Setup 

%First, check for a TAW dataset 
if strcmp(RunupMethod.Method,'TAW') || strcmp(RunupMethod.ExceedMethod,'TAW')

%% Initial parameter selection and calculations        
% Constants
g = 9.81;                       % Gravity (SI) [m/s2]
rho = 1025;                     % Density [kg/m3];
gamma = 0.78;                   % Saturated Breaking Coefficient

HS = Wavecon.Hs; %wave data
PSP = Wavecon.T; % period data
wdepth=Wavecon.z0(1);
wdepth=repmat(wdepth,length(HS),1);
if wdepth(1)>10;
    [Ho]=BackCalcHo_Limits(PSP,wdepth(1),HS);
else
    Ho=HS;
end

HS=Ho;

% Compute wave properties according to LWT   
L_STK = (g*(PSP.^2))/(2*pi);
Tm = PSP/1.1;                           % Eqn D.4.5-26

% Deep water wave length for the DIM method to be applied in the presence
% of structures/barriers.
L = (g*(Tm.^2))/(2*pi);    


Ej=SlopeInfo.toe_ele;
SWL = SlopeInfo.MSL+tide;
S=RunupMethod.StockSlope; %average regional slope for the profile

%% Begin Calculation 

%Wall slope is the maximum slope between the toe and the crest, the maximum
%slope fo a cliff wall for instance
wall_slope=(SlopeInfo.depth(SlopeInfo.overtop_point)-SlopeInfo.depth(1))/(SlopeInfo.overtop_point-1);
wall_slope=atan(wall_slope)*(180/pi);

%Get Dynamic Water level
%Bathy interp for slopes that do not have a discernable beach profile
%before the rocky outcrop/cliff/bluff 

CoastX=SlopeInfo.xpf(1);
CoastY=SlopeInfo.ypf(1);
CoastZ=SlopeInfo.depth(1);

OffshoreX=Wavecon.x0;
OffshoreY=Wavecon.y0;
OffshoreZ=Wavecon.z0-SlopeInfo.MSL; %Referenced to MSL so subtract MSL to get NAVD88
%Now calculate the distance between the points
xcomp=(CoastX-OffshoreX)^2;
ycomp=(CoastY-OffshoreY)^2;
distance=round(sqrt(xcomp+ycomp));

%Interpolate between the two 
x=[1 distance];
y=[-1*OffshoreZ CoastZ];    
newx=1:distance;
int_bath=interp1(x,y,newx);

if isempty(int_bath)
%use locally derived slope to get the new length
loc_slope=(SlopeInfo.depth(3)-SlopeInfo.depth(1))/(3-1);
x1=round((((CoastZ-(SlopeInfo.MSL))/loc_slope)-1)*-1);
x1=[x1 1];
y=[SlopeInfo.MSL CoastZ];
newx=[x1:1];
int_bath=interp1(x1,y,newx);
int_bath(end)=[];
distance=length(int_bath);

if length(distance)>0
CoastZ1=int_bath(1);
x1=round((((CoastZ1-(-1*OffshoreZ))/S)-1)*-1);
x1=[x1 1];
y=[-1*OffshoreZ CoastZ1];
newx=[x1:1];
int_bath2=interp1(x1,y,newx);
int_bath2(end)=[];
distance2=length(int_bath2);

int_bath=[int_bath2 int_bath];
distance=sum([distance distance2]); %for troubleshooting while debugging
else

x1=round((((CoastZ-(-1*OffshoreZ))/S)-1)*-1);
x1=[x1 1];
y=[-1*OffshoreZ CoastZ];
newx=[x1,1];
int_bath=interp1(x1,y,newx);
int_bath(end)=[];
distance=length(int_bath);
end
else
int_bath(end)=[];
end
%This is now the profile used to calculate the +- Hmo slope

%% Modify interpolated bathymetry
int_bath=[int_bath SlopeInfo.depth];

    %Extend the profile down to MSL 
    
    %first, get the local slope based on the first 3 points
    tmp_slope=(SlopeInfo.depth(3)-SlopeInfo.depth(1))/2;
    U=SlopeInfo.xpf(3)-SlopeInfo.xpf(1);
    V=SlopeInfo.ypf(3)-SlopeInfo.ypf(1);
    [r,az] = pcoord(U,V);
    r=100;[U,V] = rcoord(r, az);
    xpf_tmp=SlopeInfo.xpf(1)-U;
    ypf_tmp=SlopeInfo.ypf(1)-V;
    
    xpf_tmp=linspace(xpf_tmp,SlopeInfo.xpf(1),100);
    ypf_tmp=linspace(ypf_tmp,SlopeInfo.ypf(1),100);
    
    %compute distance (should actually be 100 because xshore resolution is 1m) 
    distance=100;
      
    y2=-1*abs(SlopeInfo.depth(1)+tmp_slope*(100-1));
    x=[1 distance];
    y=[y2 CoastZ]; 
    newx=1:distance;
    int_bath_tmp=interp1(x,y,newx);
    int_bath_tmp(end)=[];
    int_bath_tmp=[int_bath_tmp SlopeInfo.depth];
    int_bath_xpf=[xpf_tmp(1:end-1) SlopeInfo.xpf];
    int_bath_ypf=[ypf_tmp(1:end-1) SlopeInfo.ypf];
    
    [Tmsl,~]=curveintersect([1:length(int_bath_tmp)],int_bath_tmp,[1 length(int_bath_tmp)],[SlopeInfo.MSL SlopeInfo.MSL]);
   
    xpf_test=interp1([1:length(int_bath_tmp)],int_bath_xpf,Tmsl);
    ypf_test=interp1([1:length(int_bath_tmp)],int_bath_ypf,Tmsl);
    Tmsl=round(Tmsl(1)); 
    
        
    new_prof=int_bath_tmp(Tmsl:end);
    new_xpf=int_bath_xpf(Tmsl:end);
    new_ypf=int_bath_ypf(Tmsl:end);
    
    %Now calculate the profile the old way
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Snippet used for diagnostics while debugging
    CoastX=new_xpf(1);
    CoastY=new_ypf(1);
    CoastZ=new_prof(1);
    

    OffshoreX=Wavecon.x0;
    OffshoreY=Wavecon.y0;
    OffshoreZ=Wavecon.z0-SlopeInfo.MSL; %Referenced to MSL so subtract MSL to get NAVD88
    int_bath=[int_bath new_prof];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    S=fslope(:,3); % Set slope to a fill slope using regional average
%New int_bath calculation based on regional beach slope
   
   distance=200;
   x=[1 distance];
   OffshoreZ=S*(distance)*(-1)+CoastZ;
   y=[OffshoreZ CoastZ];
   newx=1:distance;
   int_bath=interp1(x,y,newx);
   int_bath(end)=[];
   int_bath=[int_bath new_prof];


%% Static Setup, Infragravity Swash, Incident Swash DWL and Hmo (**STK**)
%code snippet directly from Allan and others (2015)
% Apply the Stockdon et al. (2006) method for computing Static Setup
% and Infragravity Swash.
% Static Setup, Stockdon (2006) Eqn: 10
STK_setup = 0.35*S*sqrt(HS.*L_STK);    
 % Infragravity Swash, Stockdon (2006) Eqn: 12
STK_IG = 0.06*sqrt(HS.*L_STK);   
% Incident Swash, Stockdon (2006) Eqn: 11
STK_INC = 0.75*S*sqrt(HS.*L_STK);
% Stockdon (2006) Eqn: 7, 11, and 12
STK_swash = sqrt(STK_IG.^2 + STK_INC.^2)/2;
% Runup Stockdon (2006) Eqn: 9
STK_Runup = 1.1*(STK_setup + STK_swash);
% Total Water Levels
STK_TWL = STK_Runup+SWL;%Reference to NAVD88 and add in the tides
% Determine the DWL2%
STK_DWL2 = 1.1*(STK_setup + STK_IG/2) - Ej + SWL; %<----- just the DWL above the toe 

if max(STK_DWL2)<0 && strcmp(RunupMethod.Method,'TAW');
    Ej=SlopeInfo.depth(1);
    STK_DWL2 = 1.1*(STK_setup + STK_IG/2) - Ej + SWL;
end

% Wave heigth at toe calculated using a breaker index of 0.78
STK_Hmo = STK_DWL2*0.78;    
% If the depth limited breaking is larger than the offshore conditions,
% then the latter will be used.
STK_Hmo(STK_Hmo>HS) = HS(STK_Hmo>HS); % v0.1.4

%If the wave height at the toe is less than 0 then NaN out (NTC 16 Sep
%2013)
STK_Hmo(STK_Hmo<0) = NaN;

% SETUP
% v0.1.5

setup = STK_setup;
DWL2 = STK_DWL2;
Hmo = STK_Hmo;
%SWL2=SWL+(1.1*STK_setup);
SWL2=SWL;


%% Calculate local structure slope

%Round 1: calculate the slope based on the SWL +/- Hmo, so need the Hmo
%percentiles

Up=STK_TWL;
Down=SWL2;
rem=find(DWL2>0);
    
tup=Up(rem);
    
%if UpY exceeds the nearshore maxima or the clifftop elevation, set the
%elevation to be:   
Up(Up>SlopeInfo.depth(SlopeInfo.overtop_point))=SlopeInfo.depth(SlopeInfo.overtop_point)-0.01;
Down(Down>SlopeInfo.depth(SlopeInfo.overtop_point))=SlopeInfo.depth(SlopeInfo.overtop_point)-0.02;
Up2=Up;
Down2=Down;

% Set up an interative approach using linspace
Up=linspace(min(Up),max(Up),100);
Down=linspace(min(Down),max(Down),100);

Iup=NaN(length(Up),1);
Idown=NaN(length(Down),1);

UpX=[1 length(int_bath)];
DownX=[1 length(int_bath)];  
genX=[1:length(int_bath)];

%Simplify the Bathymetry first
int_bath2=int_bath;
int_bath2(int_bath2<min(Down)-1)=NaN;
oo=find(int_bath2>max(Up));
oo=oo(2:end);
int_bath2(oo)=NaN;

for ii=1:length(Up);
    UpY=[Up(ii) Up(ii)];
    DownY=[Down(ii) Down(ii)];
    if isnan(Iup(ii));
    [Tup,~]=curveintersect(genX,int_bath,UpX,UpY);
    if isempty(Tup)
        Tup=1;
    end
    Iup(ii)=Tup(1);
    [Tdown,~]=curveintersect(genX,int_bath,DownX,DownY);
    if isempty(Tdown)
        Tdown=1;
    end
    Idown(ii)=Tdown(1);
    end
end

%Bin the data for the upper and lower bounds for iteration 
%UP
if Up(1)>0;
edges=0;
edges=[edges;(edges(end)+Up(1))/2];
for ii=2:length(Up)-1;
    temp=(Up(ii)+Up(ii+1))/2;
    edges=[edges;temp];
end
edges=[edges;Up(end)];
else
   edges=Up(1)-0.01; 
   edges=[edges;(edges(end)+Up(1))/2];
   for ii=2:length(Up)-1;
    temp=(Up(ii)+Up(ii+1))/2;
    edges=[edges;temp];
   end
edges=[edges;Up(end)];
end
edgesup=edges;
clear edges
%DOWN
if Down(1)>0;
edges=0;
edges=[edges;(edges(end)+Down(1))/2];
for ii=2:length(Down)-1;
    temp=(Down(ii)+Down(ii+1))/2;
    edges=[edges;temp];
end
edges=[edges;Down(end)];
else
   edges=Down(1)-0.01; 
   edges=[edges;(edges(end)+Down(1))/2];
   for ii=2:length(Down)-1;
    temp=(Down(ii)+Down(ii+1))/2;
    edges=[edges;temp];
   end
edges=[edges;Down(end)];
end
edgesdown=edges;
clear edges

temp=Up2;
temp(temp<0)=0;
temp(temp>max(Up))=max(Up);

YUp = discretize(temp,edgesup);


temp=Down2;
temp(temp<0)=0;
temp(temp>max(Up))=max(Up);

YDown = discretize(temp,edgesdown);

NUp=Iup(YUp);
NDown=Idown(YDown);

SLOPE=abs((Up2-Down2))./abs(NUp-NDown);

%With this slope, need to now do the first iteration of the TAW runup
%code

%% TAW Reduction factors

%If assuming that this is only along plunging rocky cliffs
yr=1; %Assume a smooth rock
yb=1; %Should not have a berm in this scenario
yB=1; %No directional reduction 

%% Fix Hmo in this one instance and calculate an intial TAW runup range absed on different starting water levels
HmoNaN=find(isnan(Hmo));
Hmo_New=Hmo;Hmo_New(HmoNaN)=0.00001;

Ib_local = SLOPE./(sqrt(Hmo_New./L)); 
ybIb_local = 1.*Ib_local; %yb should only be 1 in the plunging cliff scenario
R=NaN(length(ybIb_local),1);

     temp=find(~isnan(ybIb_local));
     temp2=temp(find(ybIb_local(temp)>=0));
     temp3=temp(find(ybIb_local(temp)<1.8));
     temp4=temp(ismember(temp2,temp3));
     temp5=temp(find(ybIb_local(temp)>1.8));
     %Calculate low level and high level runup for different Ib cutoffs
     lowIb_Run=Hmo_New(temp4)*1.75*yr.*yb.*yB.*Ib_local(temp4);
     hiIb_Run=Hmo_New(temp5)*yr.*yB.*(4.3-(1.6./sqrt(Ib_local(temp5))));
     lowIb_Run=lowIb_Run + (1.1*(STK_setup(temp4)+ STK_IG(temp4)/2))+SWL(temp4);
     hiIb_Run=hiIb_Run + (1.1*(STK_setup(temp5)+ STK_IG(temp5)/2))+SWL(temp5);
     R(temp4)=lowIb_Run;
     R(temp5)=hiIb_Run;
     


%% Now to use these new runup values to define the slope  
Up=R;
Down=SWL2-(1.5*Hmo_New);

if length(Down)==length(Ej)
    Down(Down<Ej)=Ej(Down<Ej);
else
    Down(Down<Ej)=Ej;
end

    rem=find(DWL2>0);
    
    tup=Up(rem);
    
%if UpY exceeds the nearshore maxima or the clifftop elevation, set the
%elevation to be:   

Up(Up>SlopeInfo.depth(SlopeInfo.overtop_point))=SlopeInfo.depth(SlopeInfo.overtop_point)-0.01;
Down(Down>SlopeInfo.depth(SlopeInfo.overtop_point))=SlopeInfo.depth(SlopeInfo.overtop_point)-0.02;
Up2=Up;
Down2=Down;

% Set up an interative approach using linspace
Up=linspace(min(Up),max(Up),100);
Down=linspace(min(Down),max(Down),100);

Iup=NaN(length(Up),1);
Idown=NaN(length(Down),1);

UpX=[1 length(int_bath)];
DownX=[1 length(int_bath)];  
genX=[1:length(int_bath)];

%Simplify the Bathymetry first
int_bath2=int_bath;
int_bath2(int_bath2<min(Down)-1)=NaN;
oo=find(int_bath2>max(Up));
oo=oo(2:end);
int_bath2(oo)=NaN;

for ii=1:length(Up);
    UpY=[Up(ii) Up(ii)];
    DownY=[Down(ii) Down(ii)];
    if isnan(Iup(ii));
    [Tup,~]=curveintersect(genX,int_bath,UpX,UpY);
    if isempty(Tup)
        Tup=1;
    end
    Iup(ii)=Tup(1);
    [Tdown,~]=curveintersect(genX,int_bath,DownX,DownY);
    if isempty(Tdown)
        Tdown=1;
    end
    Idown(ii)=Tdown(1);
    end
end

%Bin the data for the upper and lower bounds for iteration 
%UP
if Up(1)>0;
edges=0;
edges=[edges;(edges(end)+Up(1))/2];
for ii=2:length(Up)-1;
    temp=(Up(ii)+Up(ii+1))/2;
    edges=[edges;temp];
end
edges=[edges;Up(end)];
else
   edges=Up(1)-0.01; 
   edges=[edges;(edges(end)+Up(1))/2];
   for ii=2:length(Up)-1;
    temp=(Up(ii)+Up(ii+1))/2;
    edges=[edges;temp];
   end
edges=[edges;Up(end)];
end
edgesup=edges;
clear edges
%DOWN
if Down(1)>0;
edges=0;
edges=[edges;(edges(end)+Down(1))/2];
for ii=2:length(Down)-1;
    temp=(Down(ii)+Down(ii+1))/2;
    edges=[edges;temp];
end
edges=[edges;Down(end)];
else
   edges=Down(1)-0.01; 
   edges=[edges;(edges(end)+Down(1))/2];
   for ii=2:length(Down)-1;
    temp=(Down(ii)+Down(ii+1))/2;
    edges=[edges;temp];
   end
edges=[edges;Down(end)];
end
edgesdown=edges;
clear edges

temp=Up2;
temp(temp<0)=0;
temp(temp>max(Up))=max(Up);

YUp = discretize(temp,edgesup);


temp=Down2;
temp(temp<0)=0;
temp(temp>max(Up))=max(Up);

YDown = discretize(temp,edgesdown);

NUp=Iup(YUp);
NDown=Idown(YDown);

SLOPE=abs((Up2-Down2))./abs(NUp-NDown);

SLOPE(HmoNaN)=NaN;


%% TAW Reduction factors

%yr = roughness reduciton factor that is determoined by sturcture material
%yb = berm reduction factor 
%yB = Direction reduction factor (not considered here)

yr = RunupMethod.gamma; %determined by ESI previously
if isnan(yr) 
    yr = RunupMethod.backgamma; %<--- Indicates that we are in a backshore runup equation
end

yB = 1; %Will always be 1 since we assume perpendicular waves


%%Berm width reduciton factor calculation Eqn D.4.5-21) %v0.1, v0.1.5
if SlopeInfo.berm_loc ~= 0; 

    %dh = difference of the elevated water level and the berm
    %height
    dh=STK_DWL2(rem)-SlopeInfo.depth(SlopeInfo.berm_loc);%find berm height 
 
    Bw=SlopeInfo.berm_width;

    ff=length(int_bath)-length(SlopeInfo.depth);

    tmptop = SlopeInfo.depth(SlopeInfo.berm_loc)+Hmo;
    tmpbot = SlopeInfo.depth(SlopeInfo.berm_loc)-Hmo;

    Xtop=nan(length(tmptop),1);
    Xbottom=nan(length(tmpbot),1);

for ii=1:length(tmptop);
    if isnan(tmptop(ii))
        continue
    else
        
    % Xtop (Location of STK_TWL)
            xint=find(int_bath>tmptop(ii)); xint(xint<(SlopeInfo.berm_loc+SlopeInfo.berm_width+ff))=[];xint=xint(1)-1;
            Xtop(ii) = xint;

  % Xbottom (Location of Tide + STK_setup)
            xint=find(int_bath<tmpbot(ii));xint(xint>(SlopeInfo.berm_loc+ff))=[];xint=xint(end)+1;
            Xbottom(ii) = xint;
    end
end
Lberm =  Xtop - Xbottom; %Honestly same as above
        
        
 % Determine xberm (Equation 4.5-21 & 13 (TAW report p17))JA v0.1.6
xberm=NaN(length(tmptop),1);
tmp=find((1.1*(setup(rem) + STK_swash(rem))+SlopeInfo.MSL)>-dh & -dh>0);
new_TWL=STK_TWL(rem);
tmp=find(new_TWL>-dh & -dh>0);
xberm(tmp)=new_TWL(tmp); 
tmp=find(2*Hmo(rem)>dh & dh>=0);
xberm(tmp)=2*Hmo(tmp);
xberm=fillmissing(xberm,'constant',1);
% Compute the runup reduction factor
    lb_slp=find((atan((tmptop-tmpbot)./Lberm).*(180/pi)>40)); %if this number is too high, means that the berm is too steep to count 
    yb = 1-(Bw./(2.*Lberm(rem))).*(1+cos((pi.*dh)./xberm(rem))); 
    yb(lb_slp)=1;
    gripe=NaN(length(STK_TWL),1);
    gripe(rem)=yb;
    yb=gripe;
    yb(isnan(yb)) = 1;
else
    yb=ones(length(STK_DWL2),1);
end
        
% According to the guidelines the yb must be between 0.6 and 1.
yb(yb<0.6) = 0.6;
yb(yb>1) = 1;       
          
% % Output is the composite slope and reduction factors needed to calculate the runup using TAW 
Out.Slope=SLOPE;
Out.Setup=STK_setup;
Out.yb=yb;
Out.yB=yB;
Out.yr=yr;
Out.wall_slope=wall_slope;
else 
Out=[];
   
end

end