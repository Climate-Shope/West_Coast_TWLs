function [tevent,RLs]=returnPeriods(t,data,RPs,tDclstr,mthd)
%  [tevent,RLs]=returnPeriods(t,data,RPs,tDclstr,mthd)
%  
%  INPUTS
%  t: time (datenum format)
%  data: data to be analyzed. 
%  RPs: return periods at which return values will be extracted 
%  tDclstr: time (days) between events (declustering time)
%  mthd={'LL';'GPD';'GEV'} % any combination of the three options (LL: log
%                            linear; GPD: generalized Pareto distribution; GEV: generalized extreme value)
% OUTPUTS
% RLs: return levels associated with the user defined return periods (M x N
%   matrix where the number of rows corresponds to the length of RP and the
%   number of columns corresponds to the number of methods used in the
%   analysis (i.e, each column represents the LL, GEV and or GPD in the same
%   order as input in 'mthd')
% tevent: time at which the RP occurred (limited to one time even though the same value may occur at other times) 
% note that the return period levels and dates of observations are not
% returned... should add this to the script

%% INIT
%figure('PaperPosition',[0.25 2.5 3 2.5])
hold on
rcrdlngth=(max(t)-min(t))/365; rcrdlngth=floor(rcrdlngth); %record length in years
Nspyr=1;
tevent=[]; RLs=[];
%% decluster data
[dataPks,tPks] = decluster([t data],prctile(data,50),tDclstr);
%[out] = declusterNpts([t data],tDclstr,2); tPks=out(:,1); dataPks=out(:,2);
temp=[tPks dataPks]; temp=sortrows(temp,-2);
temp=temp(1:rcrdlngth*Nspyr,:);
dataPks=temp(:,2); tPks=temp(:,1); clear temp

%% perform analyses
if sum(strcmp(mthd,'LL'))==1
    [tLL,RL_LL]=RP_LL(tPks,dataPks,RPs,rcrdlngth); 
    tevent=[tevent tLL]; RLs=[RLs RL_LL];
end

if sum(strcmp(mthd,'GPD'))==1
    [tGPD,RL_GPD]=RP_GPD(tPks,dataPks,RPs,rcrdlngth); 
    tevent=[tevent tGPD]; RLs=[RLs RL_GPD];
end

if sum(strcmp(mthd,'GEV'))==1
    [tGEV,RL_GEV]=RP_GEV(tPks,dataPks,RPs,rcrdlngth); 
    tevent=[tevent tGEV]; RLs=[RLs RL_GEV];
end

%% plot observation data onto already generated RP curves

N=length(dataPks);
rank=[1:1:N]';
P=rank./(N+1);
RpObs=1./(P*Nspyr);
semilogx(RpObs,dataPks,'ko','MarkerFaceColor',[0.5 0.5 0.5])

    
% legend
if sum(strcmp(mthd,'LL'))==1
    text(RpObs(2),max(dataPks)/3,'Log-Linear')
    plot([RpObs(2)-5;RpObs(2)-15],[max(dataPks)/3;max(dataPks)/3],'b','LineWidth',2)
end
if sum(strcmp(mthd,'GPD'))==1
    text(RpObs(2),max(dataPks)/5,'GPD')
    plot([RpObs(2)-5;RpObs(2)-15],[max(dataPks)/5;max(dataPks)/5],'g','LineWidth',2)
end
if sum(strcmp(mthd,'GEV'))==1
    text(RpObs(2),max(dataPks)/10,'GEV')
    plot([RpObs(2)-5;RpObs(2)-15],[max(dataPks)/10;max(dataPks)/10],'r','LineWidth',2)    
end
plot(RpObs(2)-5,max(dataPks)/15,'ko','MarkerFaceColor',[0.5 0.5 0.5])   
text(RpObs(2),max(dataPks)/15,'observations')
% plot(RpObs(3)-5,max(dataPks)/15,'kd')   
% text(RpObs(3),max(dataPks)/15,'return values of interest')


xlabel('Return Period (years)')

end



%% LogLinear Fit 
function [tLL,RL_LL]=RP_LL(t,data,RPfnd,rcrdlngth)
    
rank=[1:length(data)]';

rtnp=(rcrdlngth)./rank;
p=polyfit(log(rtnp),data,1);
YFIT=polyval(p,log(rtnp));

if rcrdlngth<max(RPfnd)
   RLxtrp = p(2) + p(1)*log(max(RPfnd));
   % RLxtrp = p(2) + p(1)*(max(RPfnd));

   semilogx([rtnp(1);max(RPfnd)],[YFIT(1);RLxtrp],'b--','LineWidth',1);
   hold on
end
h=semilogx(rtnp,YFIT,'b','LineWidth',2); 
grid


%  Calculate (and plot) requested return levels (RLs) using  a linear fit
for jj=1:length(RPfnd)
    RL_LL(jj,1)= p(2) + p(1)*log(RPfnd(jj,1));
    % RL_LL(jj,1)= p(2) + p(1)*(RPfnd(jj,1));

%     if RPfnd(jj)<rcrdlngth
%         delta=abs(data-RL_LL(jj,1));
%         I=find(delta==min(delta),1,'first'); tLL{jj,1}=t(I);
%     else
%         tLL{jj}=NaN;
%     end
end

semilogx(RPfnd,RL_LL,'kd','MarkerFaceColor','b','MarkerSize',8)

% find dates associated with the return period storms

for jj=1:length(RPfnd)

    if RPfnd(jj)<=rcrdlngth
        delta=abs(data-RL_LL(jj,1));
      %  tLL{jj,1}=t(find(delta>=min(delta)*0.9 & delta<=min(delta)*1.1));% find all dates within 10% of RL
        tLL{jj,1}=t(find(delta==min(delta),1,'last'));
    else
        tLL{jj}=NaN;
    end
end
end

%% GPD method
function [tGPD,RL_GPD]=RP_GPD(t,data,RPfnd,Nyrs)
    Nspyr=1; % number of samples per year for extreme value analysis

    thrsh=min(data); % threshold
    data2=data-thrsh; data2=data2(data2>0);% threshold needs to be subtracted from data prior to running GPD
    [pars,cipars]=gpfit(data2);
    k=pars(1); % shape parameter
    sigma=pars(2); %scale parameter
    klo=cipars(1,1); sigmalo=cipars(1,2);
    kup=cipars(2,1); sigmaup=cipars(2,2);
    
    
    % plot for given timeperiod
    Rp=linspace(1/Nspyr,Nyrs,Nyrs*Nspyr)';
    P=1-(1./(Rp*Nspyr));
    
    Rv=gpinv(P,k,sigma,thrsh); %Rv(1)=NaN;
    Rvlo=gpinv(P,klo,sigmalo,thrsh); Rvlo(1)=NaN;
    Rvup=gpinv(P,kup,sigmaup,thrsh); Rvup(1)=NaN;
    
    RpO=Rp; RvO=Rv;
    
    % plot extrapolated part
    NyrsXtr=100;
    Rp=linspace(1/Nspyr,NyrsXtr,NyrsXtr*Nspyr)';
    P=1-(1./(Rp*Nspyr));
    
    Rv=gpinv(P,k,sigma,thrsh); %Rv(1)=NaN;
    Rvlo=gpinv(P,klo,sigmalo,thrsh); Rvlo(1)=NaN;
    Rvup=gpinv(P,kup,sigmaup,thrsh); Rvup(1)=NaN;
    RpXtrGPD=Rp;
    RvXtrGPD=Rv;
    
    
    semilogx(RpXtrGPD,RvXtrGPD,'g--'); hold on
    % semilogx(Rp,Rvup,'g')
    % semilogx(Rp,Rvlo,'g')
    semilogx(RpO,RvO,'g','LineWidth',2)
    hold on
    % Pull out return values of interest from GPD fit. The array is 5, 10, 20, 50, 100 year return values.
    
    for jj=1:length(RPfnd)
        I=find(RpXtrGPD>=RPfnd(jj),1,'first');
        if jj==1
            RL_GPD(jj,1)=nanmean(RvXtrGPD(I));
        else
            RL_GPD(jj,1)=RvXtrGPD(I);
        end
    end
    
    semilogx(RPfnd,RL_GPD,'kd','MarkerFaceColor','g','MarkerSize',8)

    % find dates associated with the return period storms
    
    for jj=1:length(RPfnd)
        
        if RPfnd(jj)<=Nyrs
            delta=abs(data-RL_GPD(jj,1));
          %  tGPD{jj,1}=t(find(delta>=min(delta)*0.9 & delta<=min(delta)*1.1));% find all dates within 10% of RL
            tGPD{jj,1}=t(find(delta==min(delta),1,'last'));
        else
            tGPD{jj,1}=NaN;
        end
    end

    
end
    
%% GEV method
function [tGEV,RL_GEV]=RP_GEV(t,data,RPfnd,Nyrs)
  
   Nspyr=1; % number of samples per year for extreme value analysis
   [pars,cipars]=gevfit(data); k=pars(1); sigma=pars(2); mu=pars(3);
    
    Rp=linspace(1/Nspyr,Nyrs,Nyrs*Nspyr)';
    P=1-(1./(Rp*Nspyr));
    
    Rv=gevinv(P,k,sigma,mu); %Rv(1)=NaN;
    RpO=Rp; RvO=Rv;
    
    % plot original timelength and extrapolated part
    NyrsXtr=100;
    RpXtrGEV=linspace(1/Nspyr,NyrsXtr,NyrsXtr*Nspyr)';
    P=1-(1./(RpXtrGEV*Nspyr));
    RvXtrGEV=gevinv(P,k,sigma,mu); 
    semilogx(RpXtrGEV,RvXtrGEV,'r--')
    % plot up to time of real data
    semilogx(RpO,RvO,'r','LineWidth',2)
    
    % Pull out return values of interest from GEV fit. The array is 5, 10, 20, 50, 100 year return values.
    for jj=1:length(RPfnd)
        I=find(RpXtrGEV>=RPfnd(jj),1,'first');
        if jj==1
            RL_GEV(jj,1)=nanmean(RvXtrGEV(I));
        else
            RL_GEV(jj,1)=RvXtrGEV(I);
        end
    end
    semilogx(RPfnd,RL_GEV,'kd','MarkerFaceColor','r','MarkerSize',8)
        
    % find dates associated with the return period storms
    for jj=1:length(RPfnd)
        
        if RPfnd(jj)<=Nyrs
            delta=abs(data-RL_GEV(jj,1));
           % tGEV{jj,1}=t(find(delta>=min(delta)*0.9 & delta<=min(delta)*1.1));% find all dates within 10% of RL
            tGEV{jj,1}=t(find(delta==min(delta),1,'last'));
        else
            tGEV{jj,1}=NaN;
        end
    end


end



