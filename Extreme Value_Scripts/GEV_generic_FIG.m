function [res,GEV_RL,CDF_prob,CDF,GEV_conf]=GEV_generic_FIG(data,yr,return_period_array,type)
% [res,GEV_RL]=GEV_GENERIC(data,yr)
% Performs standard AMM GEV analysis on any dataset
% INPUTS: 
%   data : a vector of your annual maxima (e.g., water level, TWL, runup)
%   yr:   a vector of years for this time series
%OUTPUTS: 
%   res: all the information from the GEV analysis, where res.par_ests are
%        the GEV parameter estimates. In this case, res.par_ests(1) = xi 
%        (shape parameter), res.par_ests(2 ) = sigma (scale parameter), and
%        res.par_ests(3) = mu (location parameter).
%   GEV_RL: the GEV return levels for the 1, 10, 25, 50, 100, 500, and 1000 yr events. 
%  finalized by Ruggiero on 10/15/2007
%  Edited by K. Serafin on 3/11/2014
%  Edited by James Shope on 8/12/2019
%%

%First grab the zscores needed to calculate PDF
Probs=linspace(0,1,1000);
alpha=1-Probs;
alpha2=alpha/2;
Probs2=1-alpha2;
z = @(p) -sqrt(2) * erfcinv(p*2);
zscore = z(Probs2);


%%%%%%%%%%%%%%%%%% Annual first
which_years_ind=find(data>0);%all years, usual analysis
sortdata=sort(data(which_years_ind),'ascend'); 
sortdata_d=sort(data(which_years_ind),'descend'); 
data2=[data type];
data2=sortrows(data2,-1);
type=data2(:,2);

[res,theta]=gev(sortdata,1);

%     'Using Maximum Likelihood Method'
    k=res.par_ests(1);
    s=res.par_ests(2);
    m=res.par_ests(3);
    cov=res.varcov;


tmp(:,1)=yr;
tmp(:,2)=data;


cov2=zeros(3,3); % Reordering the covar
cov2(1,1)=cov(3,3);
cov2(1,2)=cov(1,2);
cov2(1,2)=cov(3,2);
cov2(1,3)=cov(3,1);
cov2(2,1)=cov(2,3);
cov2(2,2)=cov(2,2);
cov2(2,3)=cov(2,1);
cov2(3,1)=cov(1,3);
cov2(3,2)=cov(1,2);
cov2(3,3)=cov(1,1);

se=(diag(cov2)).^.5; %Recalc se -->already included above
ci = 1.96*se;

% Allocate array
probdata = nan(length(sortdata),1);
for ii=1:length(sortdata);
    probdata(ii)=ii/(length(sortdata)+1);  
end
ndata=-log(1-probdata);

%% estimate the confidence intervals using the delta method

count=0;

for ii=1:.1:1000;
    count=count+1;
    boink(count)=ii;
    yp(count)=-log(1-1/boink(count));
    delzpt=[1,-k^(-1)*(1-yp(count)^(-k)),s*k^(-2)*(1-yp(count)^(-k))-s*k^(-1)*yp(count)^(-k)*log(yp(count))];
    varzp=delzpt*cov2*transpose(delzpt);
    year_varzp(count)=(varzp)^.5;
    retgev=(1-(-log(1-1/boink(count)))^-k)*s/-k+m;
    retgevind(count)=retgev;
    cidelta(count)=1.96*(varzp)^.5;
end


%%%%%%%%%%%%%%%%%
%%Diagnostics
%%%%%%%%%%%%%%%%%
for ii=1:length(sortdata);
    gemp(ii)=ii/(length(sortdata)+1);
    gevprob(ii)=exp(-(1+k*(sortdata(ii)-m)/s).^-(1/k));   %Katy's way
    dum=-log(gemp(ii));
    gevneg1(ii)=m-s/k*(1-dum^-k);
    gevest(ii)=exp(-(1+k*(ii-m)/s).^-(1/k)); %Katy's way
end;


%%prob plot
figure(1)
subplot(2,2,1)
plot(gemp,gevprob,'r.')
hold on
plot(0:.1:1,0:.1:1,'k')
xlabel('empirical')
ylabel('Model')
title('probability plot')



%% quantile plot
subplot(2,2,2)
plot(gevneg1,sortdata,'r.')
hold on
plot(min(sortdata):.1:max(sortdata),min(sortdata):.1:max(sortdata),'k')
axis([gevneg1(1) gevneg1(end) gevneg1(1) gevneg1(end)])
ylabel('empirical')
xlabel('Model')
title('quantile plot')

%gev_1=((1-(-log(1-1/1))^-k)*s/-k+m)*100/100;
gev_2=((1-(-log(1-1/2))^-k)*s/-k+m)*100/100;
gev_5=((1-(-log(1-1/5))^-k)*s/-k+m)*100/100;
gev_10=((1-(-log(1-1/10))^-k)*s/-k+m)*100/100;
gev_25=((1-(-log(1-1/25))^-k)*s/-k+m)*100/100;
gev_50=((1-(-log(1-1/50))^-k)*s/-k+m)*100/100;
gev_100=((1-(-log(1-1/100))^-k)*s/-k+m)*100/100;
gev_500=((1-(-log(1-1/500))^-k)*s/-k+m)*100/100;
gev_1000=((1-(-log(1-1/1000))^-k)*s/-k+m)*100/100;

zcheck=min(data):1:gev_1000+1;
check=(1+k*(zcheck-m)./s);
bad=find(check(find(isnan(check)==0))<=0);
length(bad)


%% return period plot
subplot(2,2,3)
plot(gev_10,gev_10,'w')
hold on
plot(gev_10,gev_25,'w')
plot(gev_10,gev_50,'w')
plot(gev_10,gev_100,'w')
plot(gev_10,gev_500,'w')


plot(1./ndata,sortdata_d,'k.','MarkerSize',16)

plot(1./yp,retgevind,'r')
plot(1./yp,retgevind-cidelta,'g--')
plot(1./yp,retgevind+cidelta,'g--')

ylabel('Return Level')
xlabel('Return Period')
set(gca,'XScale','log')
title('Return Period plot')
xlim([10^-1 10^3])
%ylim([0.2 0.6])



%% pdf
%f = wgevpdf(0:.01:max(sortdata)+5,k,s,m);
% Katy's fix - looks better in quantile plots
f = gevpdf(0:.01:max(sortdata)+5,k,s,m);

gog=zeros(length(sortdata));
numbins = 8;
% Density Plot
[N,X] = hist(sortdata,numbins);
subplot(2,2,4);hold on
bar(X,N/length(sortdata)*numbins/(max(sortdata)-min(sortdata)))
plot(0:.01:max(sortdata)+5,f,'g')
plot(sortdata,gog,'k.')
title('Density plot')
ylabel('f(z)')
xlabel('z')
axis([(min(sortdata) - (max(sortdata)-min(sortdata))/10)  (max(sortdata) + (max(sortdata)-min(sortdata))/10) 0 (max(f)+.05*max(f))])
set(gcf,'Color','w')




figure(2)
%% return period plot
%subplot 211
set(gca,'FontSize',12)

plot(gev_10,gev_10,'w')
hold on
plot(gev_10,gev_25,'w')
plot(gev_10,gev_50,'w')
plot(gev_10,gev_100,'w')
plot(gev_10,gev_100,'w')
plot(gev_10,gev_100,'w')
plot(gev_10,gev_100,'w')
plot(gev_10,gev_500,'w')
plot(1./ndata,sortdata_d,'w.','MarkerSize',16)


plot(1./yp,retgevind,'r','LineWidth',3)
plot(1./yp,retgevind-cidelta,'g--','LineWidth',3)
plot(1./yp,retgevind+cidelta,'g--','LineWidth',3)
plot(1./ndata,sortdata_d,'k.','MarkerSize',16)

gg=find(type==1);
qq=1./ndata;
plot(qq(gg),sortdata_d(gg),'g.','MarkerSize',16)

gg=find(type==2);
qq=1./ndata;
plot(qq(gg),sortdata_d(gg),'b.','MarkerSize',16)

gg=find(type==3);
qq=1./ndata;
plot(qq(gg),sortdata_d(gg),'c.','MarkerSize',16)



gev_2 = gev_2*10/10;
gev_5 = gev_5*10/10;
gev_10 = gev_10*10/10;
gev_25 = gev_25*10/10;
gev_50 = gev_50*10/10;
gev_100 = gev_100*10/10;
gev_500 = gev_500*10/10;
GEV_RL = [gev_2,gev_5,gev_10,gev_25,gev_50,gev_100,gev_500, gev_1000];
%RPs=[2,5,10,25,50,100,500,1000];

%% Grab the retrun level values

for ii=1:length(return_period_array);
    clump=((1-(-log(1-1/return_period_array(ii)))^-k)*s/-k+m)*100/100;
    clump=clump*10/10;
    GEV_RL(ii)=clump;
end


%OKay search for the same values in the return plot and grab the same
%return level

for ii=1:length(return_period_array)
    rugged=GEV_RL(ii);
    [~,ind]=min(abs(rugged-retgevind));
    
    hi=rugged+zscore*year_varzp(ind);
    lo=rugged-zscore*year_varzp(ind);
    dist=[fliplr(lo),hi(2:end)];
    CDF(ii,:)=dist;
end
CDF_prob=[[1-fliplr(Probs2)],Probs2(2:end)];


[~,ia,~] = intersect(round(1./yp),round(return_period_array));

GEV_conf = cidelta(ia);

ylabel('Return Level (m)')
title('Annual Maximum Method Return Period (years)')
set(gca,'XScale','log')
xlim([0 500])
ylim([min(sortdata) max(sortdata)+0.2])



eval(['h=legend(''10-year event = ' num2str(round(gev_10*10)/10) ' m '',''25-year event = ' num2str(round(gev_25*10)/10) ' m  '',''50-year event = ' num2str(round(gev_50*10)/10) ' m '',''100-year event = ' num2str(round(gev_100*10)/10) ' m '')']);
set(h,'FontSize',10)
legend boxoff
set(h,'location','southeast')



