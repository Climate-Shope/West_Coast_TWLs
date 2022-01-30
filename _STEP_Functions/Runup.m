function [R,out1,out2] = Runup(H,T,m,method,option)
% [R,out1,out2,out3] = Runup(H,T,m,method,option) 
%
% R = runup (m above MSL)
%
% H: at toe of structure OR deepwater wave height depending on method used(m) 
%                     (Nx1 or 1x1 matrix); only needed if methods 1 are used
%
% T: deepwater period (s) (matrix size must be same as for Ho); only needed 
%                      if methods 1 are used
%
% m: foreshore slope (rads) (Nx1 or 1x1 => does not need be same size as H)
% 
% option: option for select methods:  
%   if selected method 10 (Mase) then must give Rmax, Rs, or Rmean (1,2, or 3)
%   method of calculation
% 	  1	 Larson and Kraus, 1989	** R/Ho=1.47*(m/sqrt(Ho/Lo))^0.79 **
%     2  Holman, 1986				** R=0.45*Ho+1.21   Ho = deep water wave height**
%     3  Hunt,1959&Battjes,1974     ** R=tan(beta)*sqrt(Ho/Lo)
%     4  Raubenheimer&Guza, 1996 ** R=Ho*((Iro^2)/pi) or R=Ho*sqrt(pi/(2*m)) depending if it's a reflective or saturated beach
%     5  LWT Nielsen&Hanslow, '91** non-breaking waves Not in yet... R=mWL+H((pi/2)/tanB)^0.5   
%           where H is the wave height when H=h
%     6  Walton & Ahrens, 1989	** R=(sin(m)*H)/sqrt(H/Lo);  breaking waves; H is at the foot of the structure
%     7  Walton & Ahrens, 1989	** R=H*(2*pi)^0.5*(pi/2*m)^1/4; non-breaking waves; H is at the foot of the structure
%     8  Hunt, 1959					** R=H*2.3*(tan(m)/sqrt(H/T^2)); 
%     9  Erikson, 2000				** R=0.404*((g(H*m)^0.5)/(H/L)^0.25)^2
%    10  Mase, et al., 1984         **R=Ho*d((tanB)/sqrt(Ho/Lo))^e where the coefficeints d and e are dependent on whether one is
%                                   looking for Rmax, R(1/3), or R mean
%    11  test...de Waal and van der Meer, 1992
%                                   ** R=(2ms/m+s)*sqrt(Ho/Lo), enter s as 'option' term 
%    12  Mayer and Kriebel, 1994     **R=m/2*(Xb-sqrt(Ho*Lo)*(-1+sqrt(1+(4*hb*sqrt(Ho*Lo))/(m*(Xb-sqrtHo*Lo))))
%                                     or R^2 + R*m*((hb/s)-sqrt(Ho*Lo) - hb*m*sqrt(Ho*Lo)=0
%                                     where m is the slope of foreshore above MSL, s the slope beneath MSL (from the shoreline to the break point), and hb, depth at breaking 
%                                     options should include [solvemethod Xb hb] or [solvemethod s Rmin Rmax hbmin hbmax dx] 
%                                     where solvemethod=1 for the first equation and 2 for the second equation; 
%                                     the second equation is solved by finding the closest value to 0
%                                     with ranges of R and Xb as given and step size given by dx
%                                     [R,out1,out2,out3]
%                                     out1: hbOpt
%                                     out2: MinEr
%    13 Stockdon et al, 2006       ** [R,Ir]=1.1(0.35m(Ho/Lo).^0.5 +0.5*
%                                     ((Ho/Lo(0.563*m^2 + 0.004)).^0.5 and
%                                     0.043(Ho/Lo)^0.5 for Ir<0.3


%clear RPlot hbPlot Er
Lo=1.56*T.^2;
g=9.81;

[Npts]=length(H);
out1=NaN;
out2=NaN;

if method==1
   if Npts>1
     for i=1:Npts
      R(i,1)=H(i,1)*1.47*(m(i,1)/sqrt(H(i,1)/Lo(i,1)))^0.79;
     end
   else
     R=H*1.47*(m./sqrt(H/Lo)).^0.79;
  end
elseif method==2
   for i=1:Npts
      R(i,1)=0.45*H(i,1)+1.21;
   end
elseif method==3
   for i=1:Npts
      R(i,1)=m(i,1)*sqrt(H(i,1)*Lo);
   end
elseif method==4
   for i=1:Npts
      Irc=((pi^3)/(2*m(i,1)))^(1/4);
      Iro=m(i,1)/(sqrt(H(i,1)/Lo(i,1)));
      if Iro>=Irc   % reflective beach
        R=H(i,1)*sqrt(pi/(2*m(i,1)));
      else          % saturated beach
        R=H(i,1)*((Iro^2)/pi);
     end
  end
elseif method==6
   for i=1:Npts
      R(i,1)=(sin(m(i,1))*H(i,1))/sqrt(H(i,1)/Lo(i,1));
   end
   
elseif method==7
   R=(H*(2*pi).^0.5).*(pi./(2*m)).^(1/4);
elseif method==8
   R=2.3*H.*(tan(m)./sqrt(H./T.^2));
elseif method==9
   term1=(H.*m).^1/2;
   term2=(H./Lo).^0.25;
   R=0.404*((term1./term2).^2);
elseif method==10
    if option==1
        d=2.319;
        e=0.771;
    elseif option==2
        d=1.497;
        e=0.695;
    elseif option==3
        d=1.085;
        e=0.678;
    end
    R=H*d*(m/sqrt(H/Lo))^e;
elseif method==11
    slpUndrWtr=option;
    R=((2*m*slpUndrWtr)/(m+slpUndrWtr))*sqrt(H*Lo);
elseif method==12
    if option(1,1)==1
        Xb=option(1,2);
        hb=option(1,3);
        sqrtHL=sqrt(H*Lo);
        R=m/2*(Xb-sqrtHL)*(-1+sqrt(1+(4*hb*sqrtHL)/(m*(Xb-sqrtHL)^2)));
        if R<0
           R=m/2*(Xb-sqrtHL)*(-1-sqrt(1+(4*hb*sqrt(H*Lo))/(m*(Xb-sqrt(H*Lo))^2)));
        end


    elseif option(1,1)==2
        s=option(1,2);  
        Rmin=option(1,3);
        Rmax=option(1,4);
        hbmin=option(1,5);
        hbmax=option(1,6);
        dx=option(1,7);
        
        ROpt=0;
        hbOpt=0;
        i=1;
        j=1;
        for R=Rmin:dx:Rmax
            for hb=hbmin:dx:hbmax
                Er(j,i)=abs(R^2 + R*m*((hb/s)-sqrt(H*Lo)) - hb*m*sqrt(H*Lo));
                RPlot(j,i)=R;
                hbPlot(j,i)=hb;
                j=j+1;
            end
            j=1;
            i=i+1;  
        end
        
        [MinErRow,MinErClmn]=find(Er==min(min(Er)));
        MinEr=Er(MinErRow,MinErClmn);
        ROpt=RPlot(MinErRow,MinErClmn);
        hbOpt=hbPlot(MinErRow,MinErClmn);
        
        % plot errors as a f(R and hb)
        %v=[1E-10;1E-9;1E-8;1E-7;1E-6;1E-5;1E-4;1E-3;1E-2];
        figure(100)
        clf
        c=contour(RPlot,hbPlot,Er);
        clabel(c)
        grid
        hold on
        xlabel('R')
        ylabel('hb')
        plot(ROpt,hbOpt,'rx')
        % end plot routine

        
        if MinEr<1
            R=ROpt;
        else
            R=[NaN];
        end
        

        
        
        out1=hbOpt;
        out2=MinEr;
    end

elseif method==13
        Ir=m./(sqrt(H./Lo));
        term1=0.35*m.*sqrt(H.*Lo);
        term2=0.5* (sqrt(H.*Lo.*(0.563*m.^2 + 0.004)));
        R=1.1.*(term1+term2);
        Idiss=find(Ir<0.3);
        if ~isempty(Idiss) && numel(H)==1;
            R(Idiss)=0.043*sqrt(H.*Lo);
        elseif  ~isempty(Idiss) && numel(H)>1;
           %R(Idiss)=0.043*sqrt(H.*Lo);
           R(Idiss)=0.043*sqrt(H(Idiss).*Lo(Idiss));
        end
        out1=Ir;
    
    
end

   
   

