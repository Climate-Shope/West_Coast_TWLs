function [Out] = CompositeSlope(Wavecon,test_con,SlopeInfo,RunupMethod,tide,fslope);

%script derived from code utilized by Allan and others (2015)

% code to calculate the composite slope in use for TAW, the function 
% utilizes the greatest wave condition at each transect to define this slope

%In this case, the slope will be set as the lccal slope defined by the
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
 L = (g*(PSP.^2))/(2*pi);                % Deepwater wave length (m)
 Hb = 0.39*(g^0.2)*(PSP.*HS.^2).^0.4;    % Calculate the wave breaker height (m)
 hb = Hb./gamma;                         % Calculate breaking depth (m)

L_STK = (g*(PSP.^2))/(2*pi);
% Mean wave period
Tm = PSP/1.1;                           % Eqn D.4.5-26

% Deep water wave lenght for the DIM method to be applied in the presence
% of structures/barriers.
L = (g*(Tm.^2))/(2*pi);    


%First figure out which toe is further inland and get its elevation (EJ)
if SlopeInfo.toe_loc>=SlopeInfo.toe_onshore && SlopeInfo.toe_loc<SlopeInfo.overtop_point;
    toe=SlopeInfo.toe_loc;
    if toe>0;
    if SlopeInfo.depth(toe) > 6 && SlopeInfo.depth(SlopeInfo.toe_onshore)<6;
        toe=SlopeInfo.toe_onshore;
    end
    end
else
    toe=SlopeInfo.toe_onshore;
    if SlopeInfo.depth(toe) < 6 && SlopeInfo.depth(SlopeInfo.toe_onshore)>6;
        toe=SlopeInfo.toe_loc;
    end
end


if toe>1;
   Ej=SlopeInfo.depth(toe);
else
   Ej=SlopeInfo.toe_ele; 
end

SWL = SlopeInfo.MSL+tide;
S=RunupMethod.StockSlope;  

CoastZ=SlopeInfo.depth(1);



%% Begin Calculation 
%Bathy interp for slopes that do not have a discernable beach profile
%before the rocky outcrop/cliff/bluff 

distance=0;
int_bath=SlopeInfo.depth;


%Change the Slope infodepth to generate a 1:1 slope
%Wall slope is the maximum slope between the toe and the crest, the maximum
%slope fo a cliff wall for instance
wall_slope=movingslope(SlopeInfo.depth,7);
wall_slope=max(wall_slope(toe:SlopeInfo.overtop_point));
wall_slope=atan(wall_slope)*(180/pi);

    
    %Calc the local slope based on the first 3 points
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
    
    CoastX=new_xpf(1);
    CoastY=new_ypf(1);
    CoastZ=new_prof(1);
    

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
    int_bath(end)=[];
    int_bath=[int_bath SlopeInfo.depth];
    
    
  
   S=fslope(:,3);

   %New int_bath calculation based on regional beach slope
   
   distance=200;
   x=[1 distance];
   OffshoreZ=S*(distance)*(-1)+CoastZ;
   y=[OffshoreZ CoastZ];
   newx=1:distance;
   int_bath=interp1(x,y,newx);
   int_bath(end)=[];
   int_bath=[int_bath SlopeInfo.depth];
   


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



%% Portion to correct if the toe is behind a higher elevation

%1. seach for a greater elevation point in front of toe 
if isempty(SlopeInfo.berm_toe) && SlopeInfo.toe_loc>0;

a=SlopeInfo.depth(1:toe);
b=find(a>SlopeInfo.depth(toe));

if ~isempty(b);
   prev=1;
   test1=1; 
   [j,jk]=max(a(b)); %so the DWL has to exceed this number to  count 
    jk=b(jk);
   test = 1.1*(STK_setup + STK_IG/2) - j + SWL; %just the DWL above the toe 
   mark=find(test<0);
   STK_DWL2(mark)=-1;  
end

end

%%
%Code snippet fomr Allan and others 2015
% Wave heigth at toe calculated using a breaker index of 0.78
STK_Hmo = STK_DWL2*0.78;    
% If the depth limited breaking is larger than the offshore conditions,
% then the latter will be used.
STK_Hmo(STK_Hmo>HS) = HS(STK_Hmo>HS); 
%If the wave height at the toe is less than 0 then NaN out (NTC 16 Sep
%2013)
STK_Hmo(STK_Hmo<0) = NaN;

%% SETUP
DWL2 = STK_DWL2;
Hmo = STK_Hmo;
SWL2=SWL+(1.1*STK_setup);

%% Calculate local structure slope

 % The slope will be computed between the TWL computed with stockdon and the
    % level of the static setup + SWL.

Up=STK_TWL;
Down=SWL2;
    rem=find(DWL2>0);
if ~isempty(rem);
        
    tup=Up(rem);
    
%if UpY exceeds the nearshore maxima or the clifftop elevation, set the
%elevation to be:    
Up(Up>SlopeInfo.depth(SlopeInfo.overtop_point))=SlopeInfo.depth(SlopeInfo.overtop_point)-0.01;
Down(Down>SlopeInfo.depth(SlopeInfo.overtop_point))=SlopeInfo.depth(SlopeInfo.overtop_point)-0.02;
Up2=Up;
Down2=Down;


% Set up an interative approach using linspace
Up=linspace(min(Up),max(Up),1000);
Down=linspace(min(Down),max(Down),1000);

Iup=NaN(length(Up),1);
Idown=NaN(length(Down),1);

UpX=[1 length(int_bath)];
DownX=[1 length(int_bath)];  
genX=[1:length(int_bath)];

% Simplify the Bathymetry first
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

else
   SLOPE=nan(length(Up),1); 
end


%% TAW Reduction factors

%yr = roughness reduciton factor that is determined by sturcture material
%yb = berm reduction factor 
%yB = Direction reduction factor (not considered here)

yr = RunupMethod.gamma; %determined by ESI previously
if isnan(yr) 
    yr = RunupMethod.backgamma; %<--- Indicates a necessary backshore runup equation
end

yB = 1; %Will always be 1 since we assume perpendicular waves 


%%Berm width reduction factor calculation Eqn D.4.5-21) %v0.1, v0.1.5
if  ~isempty(SlopeInfo.berm_loc); %SlopeInfo.berm_loc ~= 0 ||; 

    %dh = difference of the elevated water level and the berm
    %height  
    dh=STK_DWL2-SlopeInfo.depth(SlopeInfo.berm_loc);%find berm height 
    Bw=SlopeInfo.berm_width;
    
    %Alternative LBerm calc (based on dh+/-Hmo) - added NTC 12 Sep 2013    
     tmptop = SWL+Hmo;
     tmpbot = SWL-Hmo;
     
     tmptop(isnan(Hmo))=NaN;
     tmpbot(isnan(Hmo))=NaN;
    
    if exist('quality','var')
        g=1:length(tmptop);
        g(quality)=[];
        tmptop(g)=NaN;
        tmpbot(g)=NaN;
    end
    Xtop=nan(length(tmptop),1);
    Xbottom=nan(length(tmpbot),1);
    



for ii=1:length(tmptop);
    if isnan(tmptop(ii))
        continue
    else
        
    % Xtop (Location of STK_TWL)           
            xint=find(int_bath>tmptop(ii)); xint(xint<(SlopeInfo.toe_loc+ff))=[];xint=xint(1)-1;
            Xtop(ii) = xint;

     % Xbottom (Location of Tide + STK_setup)
            xint=find(int_bath<tmpbot(ii));xint(xint>(SlopeInfo.berm_loc+ff))=[];xint=xint(end)+1;
            Xbottom(ii) = xint;
    end
end

Lberm =  Xtop - Xbottom; %Honestly same as above
        
        
  % Determine xberm (Equation 4.5-21 & 13 (TAW report p17))JA v0.1.6
xberm=NaN(length(tmptop),1);
new_TWL=STK_TWL;
tmp=find(new_TWL>-dh & -dh>0);
xberm(tmp)=new_TWL(tmp);
tmp=find(2*Hmo>dh & dh>=0);
xberm(tmp)=2*Hmo(tmp);
xberm=fillmissing(xberm,'constant',1);      
  % Compute the runup reduction factor
    lb_slp=find((atan((tmptop-tmpbot)./Lberm).*(180/pi)>40)); %if this number is too high, means that the berm is too steep to count 
    yb = 1-(Bw./(2.*Lberm)).*(1+cos((pi.*dh)./xberm)); 
    yb(lb_slp)=1;
    yb(isnan(yb)) = 1;   
else
    yb=ones(length(STK_DWL2),1);  
end
        
% According to the guidelines the yb must be between 0.6 and 1.
yb(yb<0.6) = 0.6;
yb(yb>1) = 1;       
          
SLOPE(isnan(Hmo))=NaN;
yb(isnan(Hmo))=NaN;



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