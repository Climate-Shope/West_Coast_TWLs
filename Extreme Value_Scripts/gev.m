function [res,theta]=gev(data,blocksz),
%Fits generalized extreme value distribution (GEV) to block maxima data
%
%USAGE: res=gev(data,blocksz)
%
%   data: Data vector
%blocksz: Blocksize
%
%    res: Fitted distribution
%
%         res.par_ests: Estimated parameters. 1X3 vector: 
%                                                        1st element: xi
%                                                        2nd element: sigma
%                                                        3rd element: mu
%           res.funval: Value of the negative log likelihood
%       res.terminated: Termination condition. 1 if successfully terminated
%          res.details: Details of the nonlinear minimization process of the negative
%                       likelihood
%           res.varcov: Variance-covariance matrix of the parameters
%          res.par_ses: Standard deviations of the parameters of the distribution
%             res.data: Extremes of the blocks

warning off
if (~isnan(blocksz))
    n_all=length(data);
    data=block(data,blocksz,'max');
end
data=data(~isnan(data));
n=length(data);
sigma0=sqrt((6*var(data))/pi);
mu0=mean(data)-0.57722*sigma0;
xi0=0.1;
theta=[xi0,sigma0,mu0];
tmp=data;
opts=optimset('MaxFunEvals',10000,'MaxIter',2000,'TolX',1e-6,'TolFun',1e-6,'Display','off');

y = 1 + (theta(1) * (data - theta(3)))./theta(2);
	if((theta(2) < 0) | (min(y) < 0))
		theta(2)=1;
        theta(3)=1;
    end
        
[res.par_ests,res.funval,res.terminated,res.details]=fminsearch('negloglikgev',theta,opts,tmp);
[res.par_ests,res.funval,res.terminated,res.details]=fminunc('negloglikgev',res.par_ests,opts,tmp);

res.varcov = hessigev('negloglikgev',res.par_ests,data);
res.par_ses=sqrt(diag(res.varcov));
res.data=data;

warning on


        
