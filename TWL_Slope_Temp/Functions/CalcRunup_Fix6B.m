function [R,R_mag,type,DWL] = CalcRunup_Fix6B(Wavecon, SlopeInfo, RunupMethod, Composite_Slope,tide);
%Calculate Runup using slope, water level,and runup method inputs with
%slight modifications for 6B profiles

%% Initial parameter selection and calculations         (v0.1, v0.1.1)
% Constants
g = 9.81;                       % Gravity (SI) [m/s2]
rho = 1025;                     % Density [kg/m3];
gamma = 0.78;                   % Saturated Breaking Coefficient
ESI=SlopeInfo.ESI;
if isempty(ESI);
    ESI='none';
elseif isnan(ESI);
    ESI='none';
end

%Determine which toe estimation is further inland (6 m limit determiend
%through testing for the region)
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

%EJ is the elevation of the selected toe
if toe>1;
   Ej=SlopeInfo.depth(toe);
else
   Ej=SlopeInfo.toe_ele; 
end

%Still water level is ambient MSL and Tides (tides + NTRs)
SWL = SlopeInfo.MSL+tide;
%predetermined slope to use with stockdon runup method
S=RunupMethod.StockSlope; 

%If empty, means no TAW computation, but filling in the variables to ease
%computation below
if ~isempty(Composite_Slope)
yr=Composite_Slope.yr;
yr(isnan(yr))=1;
yb=Composite_Slope.yb;
yb(isnan(yb))=1;
yB=Composite_Slope.yB;
end
%% Get deepwater wave conditions to calculate the Stockdon conditions

Hs=Wavecon.Hs;
Tp=Wavecon.T;
Tp(Tp<0.1)=0.1;
wdepth=Wavecon.z0(1);
wdepth=repmat(wdepth,length(Hs),1);
%convert form shallow water/transitional conditions to deep water
if wdepth(1)>10;
    [Ho]=BackCalcHo_Limits(Tp,wdepth(1),Hs);
else
    Ho=Hs;
end

%% Dervative Wave Calculations based on Stockdon and Linear Theory
%Code snippet from Allan and others 2015
 L = (g*(Tp.^2))/(2*pi);                % Deepwater wave length (m)
 Hb = 0.39*(g^0.2)*(Tp.*Hs.^2).^0.4;    % Calculate the wave breaker height (m)
 hb = Hb./gamma;                         % Calculate breaking depth (m)

 L_STK = (g*(Tp.^2))/(2*pi);
 Tm = Tp/1.1;                           % Eqn D.4.5-26

 % Deep water wave length for the DIM method to be applied in the presence
 % of structures/barriers.
 L = (g*(Tm.^2))/(2*pi);    

 % Static Setup, Stockdon (2006) Eqn: 10
STK_setup = 0.35*S*sqrt(Ho.*L_STK);    
 % Infragravity Swash, Stockdon (2006) Eqn: 12
STK_IG = 0.06*sqrt(Ho.*L_STK);   
STK_DWL2 = 1.1*(STK_setup + STK_IG/2) - Ej + SWL; %<---- Dynamic water level over the Toe elevation
% Wave height at toe calculated using a breaker index of 0.78
STK_Hmo = STK_DWL2*0.78;    
% If the depth limited breaking is larger than the offshore conditions,
% then the latter will be used.
STK_Hmo(STK_Hmo>Ho) = Ho(STK_Hmo>Ho); % v0.1.4

%If the wave height at the toe is less than 0 then NaN out (NTC 16 Sep
%2013)
type=zeros(length(STK_Hmo),1);
type(STK_Hmo>0)=1;
STK_Hmo(STK_Hmo<0) = NaN;
Hmo=STK_Hmo;
DWL = 1.1*(STK_setup + STK_IG/2)+ SWL; 
SWL2=SWL;
if numel(Ej)>1;
    rt=Ej-SWL2;
    rt=find(rt>0);
    SWL2(rt)=Ej(rt);
else
    SWL2(SWL2<Ej)=Ej;
end



%% determine inital runup conditions
    %Use Stockdon Runup formulations 
    [R,out1,out2] = Runup(Ho,Tp,S,13); R_mag=R; R=R+SWL; %Reference it to NAVD88 by adding the Still Water level
    R_guess=R;

%% compute the exceed methodology, which should only apply where stockdon method changes to TAW
if strcmp(RunupMethod.ExceedMethod,'TAW');
    %First,determine if the DWL exceeds the toe and then if 
    %Runup exceeds the toe
    temp=find(STK_DWL2>0); %where the Taw exceedence is utilized 
    temp2=find(R>Ej);temp3=temp2(~ismember(temp2,temp));
 
    %Set up TAW calculation parameters
       
    Ib_local = Composite_Slope.Slope./(sqrt(Hmo./L)); 
    ybIb_local = yb.*Ib_local;
    
    temp=find(~isnan(ybIb_local)); 
    
     temp2=temp(find(ybIb_local(temp)>=0));
     temp3=temp(find(ybIb_local(temp)<1.8));
     
     temp4=temp(ismember(temp2,temp3));
     
     temp2a=temp(find(ybIb_local(temp)>1.8));
     temp3a=temp(find(ybIb_local(temp)<8));
     temp5=temp(find(ismember(temp2a,temp3a)));
     
     %basically if the utilized tow is <0, the use the precomputed
     %elevation
     if toe>0
     toe_d=SlopeInfo.depth(toe);
     else
       toe_d=SlopeInfo.toe_ele;  
     end
     
     
     %parameters of these conditions was determined through testing for the
     %region    
     
     if toe>=2 && ~strcmp(ESI,'6B') && Composite_Slope.wall_slope<45 && toe_d>SlopeInfo.MHW;%Composite_Slope.wall_slope>45 && toe>2
     temp6=temp((ybIb_local(temp)>8) & Ho(temp)>Hmo(temp));
     temp7=temp((ybIb_local(temp)>8) & Ho(temp)<Hmo(temp));
     
     elseif toe>=2 && ~strcmp(ESI,'6B') && Composite_Slope.wall_slope>45;
          temp6=temp(DWL(temp)>SlopeInfo.toe_ele & Ho(temp)>=Hmo(temp));
          temp7=temp(DWL(temp)>SlopeInfo.toe_ele & Ho(temp)<Hmo(temp));
          temp4=[];
          temp5=[];
     elseif Composite_Slope.wall_slope>45 && toe==-1 && ~strcmp(ESI,'6B') || toe_d<=SlopeInfo.MHW && ~strcmp(ESI,'6B');
         temp6=temp(DWL(temp)>SlopeInfo.toe_ele & Ho(temp)>=Hmo(temp));
          temp7=temp(DWL(temp)>SlopeInfo.toe_ele & Ho(temp)<Hmo(temp));
          temp4=[];
          temp5=[];
     elseif strcmp(ESI,'6B') && Composite_Slope.wall_slope>45;
          temp6=temp(DWL(temp)>SlopeInfo.toe_ele & Ho(temp)>=Hmo(temp));
          temp7=temp(DWL(temp)>SlopeInfo.toe_ele & Ho(temp)<Hmo(temp));
          temp4=[];
          temp5=[]; 
     elseif strcmp(ESI,'6B') && Composite_Slope.wall_slope<45;

    temp6=temp((ybIb_local(temp)>8) & Ho(temp)>Hmo(temp));
     temp7=temp((ybIb_local(temp)>8) & Ho(temp)<Hmo(temp));

      temp2a=temp(find(ybIb_local(temp)>1.8));
     temp3a=temp(find(ybIb_local(temp)<8));
     temp5=temp(find(ismember(temp2a,temp3a)));
     
     elseif Composite_Slope.wall_slope>45 && toe<2
         temp6=temp(SWL(temp)>SlopeInfo.toe_ele & Ho(temp)>=Hmo(temp));
          temp7=temp(SWL(temp)>SlopeInfo.toe_ele & Ho(temp)<Hmo(temp));    
     else
     temp6=temp((ybIb_local(temp)>8) & Ho(temp)>Hmo(temp));
     temp7=temp((ybIb_local(temp)>8) & Ho(temp)<Hmo(temp));
     temp2a=temp(find(ybIb_local(temp)>1.8));
     temp3a=temp(find(ybIb_local(temp)<8));
     temp5=temp(find(ismember(temp2a,temp3a)));
     end
     
     %TAW runup parameters for threhsolds defined in TAW guidance
     lowIb_Run=Hmo(temp4)*1.75*yr.*yb(temp4).*yB.*Ib_local(temp4);
     hiIb_Run=Hmo(temp5)*yr.*yB.*(4.3-(1.6./sqrt(Ib_local(temp5))));

     surge_Run=Hmo(temp6)*1.5; %SPM 1984 for runup on vertical wall 
     Walton_Run=Hmo(temp7).*sqrt(2*pi).*(pi./(2*Composite_Slope.Slope(temp7))).^(1/4); %alternative runup equation for comparison, but never was applied in practice
     
     R_DWL=R;%RDWL here means runup with DWL included, not just the runup magnitude
     %If on top of the DWL
     %R magnitude is the calculated runup magnitude above SWL, not the
     %absolute elevation relative to NAVD88
     R_mag(temp4)=lowIb_Run;
     R_mag(temp5)=hiIb_Run;
     R_mag(temp6)=surge_Run;
     R_mag(temp7)=Walton_Run;
     lowIb_Run=lowIb_Run + SWL(temp4)+(1.1*(STK_setup(temp4)));
     hiIb_Run=hiIb_Run +SWL(temp5)+(1.1*(STK_setup(temp5)));
     surge_Run=surge_Run +SWL(temp6);
     Walton_Run=Walton_Run +SWL(temp7);
     type(temp6)=2;
     type(temp7)=2;
     
     %separate parts for each non-sotckdon method runup parameters
      gg=find(lowIb_Run<DWL(temp4));
      hh=temp4(gg);
     
      lowIb_Run(gg)=R(hh); 
      type(hh)=0;
     
     gg=find(hiIb_Run<DWL(temp5));
     hh=temp5(gg);
     type(hh)=0;
      hiIb_Run(gg)=R(hh); 
      
     gg=find(surge_Run<DWL(temp6));
     hh=temp6(gg);
     type(hh)=3;
     surge_Run(gg)=R(hh);   
     
     gg=find(Walton_Run<DWL(temp7));
     hh=temp7(gg);
     type(hh)=3;
     Walton_Run(gg)=R(hh);   
      
     %classify the Runup by the potential different runup calculations
     R_DWL(temp4)=lowIb_Run;
     R_DWL(temp5)=hiIb_Run;
     R_DWL(temp6)=surge_Run;
     R_DWL(temp7)=Walton_Run;

 %The Walton runup formulation was investigated, and while incldued in the
 %output, the conditions for its applicaiton were never achieved in
 %practice and were therefore not included in the method documentation.      
    if isempty(Composite_Slope)
        type(type==1)=0;
    end
        
      R=R_DWL;
    
end

end







