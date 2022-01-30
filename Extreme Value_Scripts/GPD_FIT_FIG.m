function [GPD_RL,GPD_conf,CDF,CDF_prob,Num_Points,Percentile,res]=GPD_FIT_FIG(data,raw_length,return_period_array,type); 
%Derived from Coles 2001
%original Code from J. Allan and K. Serafin. 
%Modified by J. Shope 2019
Nums=[61,76,92,123];% Iterative number of points to include in the GPD
Count=0;


for jj=Nums;
    
    Count=Count+1;
    try
    goo_u=jj/raw_length;
    T=raw_length/365.25;
    num_exceed=jj;
    num_per_year=num_exceed/T;
    
    sortdata=sort(data);
    sortdata=sortdata(end-(jj-1):end);
    thresh=min(sortdata);
    exceedances=sortdata;
    excess=exceedances-thresh;
    xbar=mean(excess); % mean of excess
    s2=var(excess); % standard deviation of excess
    xi0=-0.5*(((xbar^2)/s2)-1); % Correction by Andrea Colombo May 9, 2005
    sigma=exp(0.5*xbar*(((xbar^2)/s2)+1));
    nu=num_per_year-1;
    theta=[nu,xi0,sigma];
    opts=optimset('MaxFunEvals',5000,'MaxIter',1000,'TolX',1e-6,'TolFun',1e-6,'Display','off');
    
    [res.par_ests,res.funval,res.terminated,res.details]=fminsearch('negloglikgpd2',theta,opts,excess,T);
    theta2=[res.par_ests(:,2) res.par_ests(3)];
    res.varcov = hessigev('negloglikgpd',theta2,excess);
    res.par_ses=sqrt(diag(res.varcov));
    res.data=sortdata;
    sigma=res.par_ests(3);
    xi=res.par_ests(2);
    nu=res.par_ests(1);
    V=zeros(3,3);
    covmatrix=res.varcov;
    
    cov2=zeros(2,2);
    cov2(1,1)=covmatrix(2,2);
    cov2(1,2)=covmatrix(2,1);
    cov2(2,1)=covmatrix(1,2);
    cov2(2,2)=covmatrix(1,1);


    V(1,1)=goo_u*(1-goo_u)/raw_length;
    V(2,2)=cov2(1,1);
    V(2,3)=cov2(1,2);
    V(3,2)=cov2(2,1);
    V(3,3)=cov2(2,2);
    
    for ii=1:length(sortdata)
        probdata(ii)=ii/(length(sortdata)+1);  
    end
    ndata=1./(1-probdata);
    
    for ii=1:length(sortdata);
    gemp(ii)=ii/(length(sortdata)+1);
    gpdprob(ii)=1-(1+xi*excess(ii)/sigma)^(-1/xi);
    quant1(ii) = ((gemp(ii)^-xi)-1);
    gpdneg1(ii)= thresh+(sigma/xi)*(quant1(ii));
  
    end
    
    
    n=.1:.1:1000;
    for ii = 1:length(n);
        xm(ii)=thresh+sigma/xi*((365.25*n(ii)*goo_u)^xi-1);
        delxT=[sigma*(365.25*n(ii))^xi*goo_u^(xi-1),...
        xi^(-1)*((365.25*n(ii)*goo_u)^xi-1),...
        -sigma*xi^(-2)*((365.25*n(ii)*goo_u)^xi-1)+sigma*xi^(-1)*(365.25*n(ii)*goo_u)^xi*log(365.25*n(ii)*goo_u)];
        varxm(ii)=delxT*V*transpose(delxT);
    
    end
    xm=xm.*100./100;
    gev_100=(thresh+sigma/xi*((365.25*100*goo_u)^xi-1))*100/100;
    gg=find(xm==gev_100);
    
    
    if ~isempty(find(varxm<0))
       CI(Count)=NaN;
    else
       CI(Count)=1.96*varxm(gg).^.5;  
    end
        
    catch
    CI(Count)=NaN;
    end
end
clear varxm xm gpdneg1 quant1 gpdprob gemp probdata
CI(isnan(CI))=999;

a=find(CI<900);
if isempty(a)
Nums=[180,185,216,247,278,309,340,371,500]; %secondary iterative approach to the number of points used in the GPD analysis
Count=0;


for jj=Nums;
    
    Count=Count+1;
    try
    goo_u=jj/raw_length;
    T=raw_length/365.25;
    num_exceed=jj;
    num_per_year=num_exceed/T;
    
    sortdata=sort(data);
    sortdata=sortdata(end-(jj-1):end);
    thresh=min(sortdata);
    exceedances=sortdata;
    excess=exceedances-thresh;
    xbar=mean(excess); % mean of excess
    s2=var(excess); % standard deviation of exvess
    xi0=-0.5*(((xbar^2)/s2)-1); % Correction by Andrea Colombo May 9, 2005
    sigma=exp(0.5*xbar*(((xbar^2)/s2)+1));
    nu=num_per_year-1;
    theta=[nu,xi0,sigma];
    opts=optimset('MaxFunEvals',5000,'MaxIter',1000,'TolX',1e-6,'TolFun',1e-6,'Display','off');
    
    [res.par_ests,res.funval,res.terminated,res.details]=fminsearch('negloglikgpd2',theta,opts,excess,T);
    theta2=[res.par_ests(:,2) res.par_ests(3)];
    res.varcov = hessigev('negloglikgpd',theta2,excess);
    res.par_ses=sqrt(diag(res.varcov));
    res.data=sortdata;
    sigma=res.par_ests(3);
    xi=res.par_ests(2);
    nu=res.par_ests(1);
    V=zeros(3,3);
    covmatrix=res.varcov;
    
    cov2=zeros(2,2);
    cov2(1,1)=covmatrix(2,2);
    cov2(1,2)=covmatrix(2,1);
    cov2(2,1)=covmatrix(1,2);
    cov2(2,2)=covmatrix(1,1);


    V(1,1)=goo_u*(1-goo_u)/raw_length;
    V(2,2)=cov2(1,1);
    V(2,3)=cov2(1,2);
    V(3,2)=cov2(2,1);
    V(3,3)=cov2(2,2);
    
    for ii=1:length(sortdata)
        probdata(ii)=ii/(length(sortdata)+1);  
    end
    ndata=1./(1-probdata);
    
    for ii=1:length(sortdata);
    gemp(ii)=ii/(length(sortdata)+1);
    gpdprob(ii)=1-(1+xi*excess(ii)/sigma)^(-1/xi);
    quant1(ii) = ((gemp(ii)^-xi)-1);
    gpdneg1(ii)= thresh+(sigma/xi)*(quant1(ii));
  
    end
    
    
    n=.1:.1:1000;
    for ii = 1:length(n);
        xm(ii)=thresh+sigma/xi*((365.25*n(ii)*goo_u)^xi-1);
        delxT=[sigma*(365.25*n(ii))^xi*goo_u^(xi-1),...
        xi^(-1)*((365.25*n(ii)*goo_u)^xi-1),...
        -sigma*xi^(-2)*((365.25*n(ii)*goo_u)^xi-1)+sigma*xi^(-1)*(365.25*n(ii)*goo_u)^xi*log(365.25*n(ii)*goo_u)];
        varxm(ii)=delxT*V*transpose(delxT);
    
    end
    xm=xm.*100./100;
    gev_100=(thresh+sigma/xi*((365.25*100*goo_u)^xi-1))*100/100;
    gg=find(xm==gev_100);
    CI(Count)=1.96*varxm(gg).^.5;
    
    catch
    CI(Count)=NaN;
    end
end


end
clear varxm xm gpdneg1 quant1 gpdprob gemp probdata
CI(isnan(CI))=999;

a=find(CI<900);
[~,ind]=min(CI);

%% Now make plots using the selected value
jj=Nums(ind);
Percentile=(1-(jj/raw_length))*100;

 goo_u=jj/raw_length;
    T=raw_length/365.25;
    num_exceed=jj;
    num_per_year=num_exceed/T;
    data2=[data type];
    data2=sortrows(data2,1);
    type=data2(end-(jj-1):end,2);
    sortdata=sort(data);
    sortdata=sortdata(end-(jj-1):end);
    thresh=min(sortdata);
    exceedances=sortdata;
    excess=exceedances-thresh;
    xbar=mean(excess); % mean of excess
    s2=var(excess); % standard deviation of exvess
    xi0=-0.5*(((xbar^2)/s2)-1); % Correction by Andrea Colombo May 9, 2005
    sigma=exp(0.5*xbar*(((xbar^2)/s2)+1));
    nu=num_per_year-1;
    theta=[nu,xi0,sigma];
    opts=optimset('MaxFunEvals',5000,'MaxIter',1000,'TolX',1e-6,'TolFun',1e-6,'Display','off');
    
    [res.par_ests,res.funval,res.terminated,res.details]=fminsearch('negloglikgpd2',theta,opts,excess,T);
    theta2=[res.par_ests(:,2) res.par_ests(3)];
    res.varcov = hessigev('negloglikgpd',theta2,excess);
    res.par_ses=sqrt(diag(res.varcov));
    res.data=sortdata;
    sigma=res.par_ests(3);
    xi=res.par_ests(2);
    nu=res.par_ests(1);
    V=zeros(3,3);
    covmatrix=res.varcov;
    
    cov2=zeros(2,2);
    cov2(1,1)=covmatrix(2,2);
    cov2(1,2)=covmatrix(2,1);
    cov2(2,1)=covmatrix(1,2);
    cov2(2,2)=covmatrix(1,1);


    V(1,1)=goo_u*(1-goo_u)/raw_length;
    V(2,2)=cov2(1,1);
    V(2,3)=cov2(1,2);
    V(3,2)=cov2(2,1);
    V(3,3)=cov2(2,2);
    
    for ii=1:length(sortdata)
        probdata(ii)=ii/(length(sortdata)+1);  
    end
    ndata=1./(1-probdata);
    
    for ii=1:length(sortdata);
    gemp(ii)=ii/(length(sortdata)+1);
    gpdprob(ii)=1-(1+xi*excess(ii)/sigma)^(-1/xi);
    quant1(ii) = ((gemp(ii)^-xi)-1);
    gpdneg1(ii)= thresh+(sigma/xi)*(quant1(ii));
  
    end
    
    
    n=.1:.1:1000;
    for ii = 1:length(n);
        xm(ii)=thresh+sigma/xi*((365.25*n(ii)*goo_u)^xi-1);
        delxT=[sigma*(365.25*n(ii))^xi*goo_u^(xi-1),...
        xi^(-1)*((365.25*n(ii)*goo_u)^xi-1),...
        -sigma*xi^(-2)*((365.25*n(ii)*goo_u)^xi-1)+sigma*xi^(-1)*(365.25*n(ii)*goo_u)^xi*log(365.25*n(ii)*goo_u)];
        varxm(ii)=delxT*V*transpose(delxT);
    
    end

%%Now for Plots
figure
%%prob plot
subplot(2,1,2)
subplot(2,2,1)
plot(gemp,gpdprob,'r.')
hold;
plot(0:.1:1,0:.1:1,'k')
xlabel('empirical')
ylabel('Model')
title('probability plot')

bigval=max([max(gpdneg1) max(sortdata)]);
smallval=min([min(gpdneg1) min(sortdata)]);

%% quantile plot

subplot(2,2,2)
plot(fliplr(gpdneg1),thresh+excess,'r.')

hold;
plot(smallval:.1:bigval,smallval:.1:bigval,'k')
ylabel('empirical')
xlabel('Model')
title('quantile plot')
axis([smallval bigval smallval bigval])


%% Return Period Plot
gev_1=(thresh+sigma/xi*((365.25*1*goo_u)^xi-1))*100/100;
gev_10=(thresh+sigma/xi*((365.25*10*goo_u)^xi-1))*100/100;
gev_25=(thresh+sigma/xi*((365.25*25*goo_u)^xi-1))*100/100;
gev_50=(thresh+sigma/xi*((365.25*50*goo_u)^xi-1))*100/100;
gev_100=(thresh+sigma/xi*((365.25*100*goo_u)^xi-1))*100/100;
gev_500=(thresh+sigma/xi*((365.25*500*goo_u)^xi-1))*100/100;
gev_1000=(thresh+sigma/xi*((365.25*1000*goo_u)^xi-1))*100/100;

subplot(2,2,3)
plot(gev_10,gev_10,'w')
hold
plot(gev_10,gev_25,'w')
plot(gev_10,gev_50,'w')
plot(gev_10,gev_100,'w')
hold on
plot(gev_10,gev_500,'w')
plot(ndata/(365*goo_u),sortdata,'k.','MarkerSize',16)

plot(n,xm,'r')
plot(n,xm-1.96*varxm.^.5,'g--')
plot(n,xm+1.96*varxm.^.5,'g--')
ylabel('Return Level')
xlabel('Return Period')
set(gca,'XScale','log')
title(['Return Period plot: Point Num ' num2str(jj)])
axis([.1 1000 thresh max(xm+nanmax(1.96*varxm.^.5))])
if max(sortdata)>gev_100;
    text(100,max(sortdata)-1,'ugh!')
else
end




%% pdf
f = gppdf(min(sortdata):.1:max(sortdata),xi,sigma,thresh);
% 
numbins=10;
% % Density Plot
 [N,X] = hist(sortdata,numbins);
 subplot(2,2,4);hold on
 bar(X,N/length(sortdata)*numbins/(max(sortdata)-min(sortdata)))
 plot(min(sortdata):.1:max(sortdata),f,'g')
 title('Density plot')
 ylabel('f(z)')
 xlabel('z')

 set(gcf,'Color','w');
 
 Out='done';
 
 
 
 %% separate Retun Period Plot
 
 figure;
plot(gev_10,gev_10,'w')
hold
plot(gev_10,gev_25,'w')
plot(gev_10,gev_50,'w')
plot(gev_10,gev_100,'w')
hold on
plot(gev_10,gev_500,'w')
plot(ndata/(365*goo_u),sortdata,'k.','MarkerSize',16)

gg=find(type==1);
qq=ndata/(365*goo_u);
plot(qq(gg),sortdata(gg),'g.','MarkerSize',16)
gg=find(type==2);
qq=ndata/(365*goo_u);
plot(qq(gg),sortdata(gg),'b.','MarkerSize',16)
gg=find(type==2);
qq=ndata/(365*goo_u);
plot(qq(gg),sortdata(gg),'c.','MarkerSize',16)
plot(n,xm,'r')
plot(n,xm-1.96*varxm.^.5,'g--')
plot(n,xm+1.96*varxm.^.5,'g--')
ylabel('Return Level')
xlabel('Return Period')
set(gca,'XScale','log')
title(['Return Period plot:Point Num ' num2str(jj)])
axis([.1 1000 thresh max(xm+nanmax(1.96*varxm.^.5))])
if max(sortdata)>gev_100;
    text(100,max(sortdata)-1,'ugh!')
else
end 
 
%%Spit out CDF for Confidence Intervals and Probability Topping then
 %ultimately save the one with the better constraints 
 n=round(n,1);
 ind=find(ismember(n,round(return_period_array)));
 GPD_conf=1.96*varxm(ind).^.5;
 GPD_RL=(thresh+sigma/xi*((365.25*return_period_array*goo_u).^xi-1))*100/100;
 GPD_RL=GPD_RL';
 
 %Now construct the CDF for this one
 %First grab the zscores needed to calculate PDF
    Probs=linspace(0,1,1000);
    alpha=1-Probs;
    alpha2=alpha/2;
    Probs2=1-alpha2;
    z = @(p) -sqrt(2) * erfcinv(p*2);
    zscore = z(Probs2);
 
 for ii=1:length(return_period_array)
    rugged=GPD_RL(ii);
    send=ind(ii);
    hi=xm(send)+zscore.*varxm(send).^.5;
    lo=xm(send)-zscore.*varxm(send).^.5;
    dist=[fliplr(lo),hi(2:end)];
    CDF(ii,:)=dist;
 end
CDF_prob=[[1-fliplr(Probs2)],Probs2(2:end)];
Num_Points=jj;
