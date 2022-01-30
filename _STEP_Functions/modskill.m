function skill=modskill(pred,obs,dim)
% MSKILL - Calculate model performance or "skill" measures
%
%   SKILL =  MSKILL(PRED,OBS) calculates model performance measures
%   suggested by Willmott (1982) betweeen predicted or modeled 
%   parameter, PRED, and observed value, OBS. NaNs are removed. 
%   PRED and OBS should have the same size and shape. If a matrix 
%   is provided for PRED and OBS, then the model skill will be 
%   computed for each column-wise (unless an optional DIM argument 
%   is provided).
% 
%   SKILL = MSKILL(...,DIM) - computes the model skill along the 
%   dimension DIM.  MSKILL is only valid for 2D matrices, so DIM 
%   should either be 1 (default) or 2.
%
%   OUTPUT - Function returns a (possibly n dimensional structure):
%       
%       'n'         -  number of elements 
%       'meanp'     -  mean of predicted values
%       'stdp'      -  standard deviation of predicted values
%       'meano'     -  mean of observations
%       'stdo'      -  standard deviation of predicted values
%       'slope'     -  slope of least squares linear line
%       'intercept' -  intercept of least squares linear line 
%       'rmse'      -  root mean square error added by LHE
%       'rmsa'      -  systematic error
%       'rmsu'      -  unsystematic error
%       'd'         -  index of agreement
%
% REFERENCE 
%   Willmot, C.J., (1982) Some comments on the evaluation
%      of model performance, Bulletin of American Meteorological 
%      Society, vol. 63, pp 1309-1313.
%
% EXAMPLE
% %Create some data
%  x =0:0.1:2*pi;
%  pred = sin(x);
%  obs = pred+0.2*rand(1,numel(pred));
%  
% %calculate the skill of the estimate
%  skill = mskill(pred,obs);


%check inputs
if (~isequal(numel(pred),numel(obs)) || ...
        ~isequal(size(pred),size(obs)));
    error('Inputs should be vectors or matrices of the same size');
end
if ~exist('dim','var')
   dim=1; 
else
    if dim==0 || dim>2 
        error('Dimension can only be 1 or 2');
    end
end
[m,n,z]=size(obs);
if z>1
    error('MSKILL only works on 2-D matrices');
end

if m==1 || n==1 %vector case
    skill = calc_error(pred,obs);
else %matrix case
    predc=num2cell(pred,dim);
    obsc=num2cell(obs,dim);
    skill=cellfun(@(x,y)(calc_error(x,y)),predc,obsc);
end

%engine-------------------------------------------------------------
function skill=calc_error(pred,obs)
%remove nans
idx = isnan(pred) | isnan(obs);
pred(idx)=[];
obs(idx)=[];

%summary measures
skill.n=numel(pred); 
skill.meanp=mean(pred);
skill.stdp=std(pred);
skill.meano=mean(obs);
skill.stdo=std(obs);

%linear correlation
p=polyfit(obs,pred,1);
skill.slope=p(1);
skill.intercept=p(2);

%difference measures
phat=polyval(p,obs);
skill.rmse=sqrt((sum((pred-obs).^2)/skill.n)); %root-mean square error
skill.rmsa=sqrt((sum((phat-obs).^2))./skill.n); %systematic error
skill.rmsu=sqrt((sum((pred-phat).^2))./skill.n); %unsystematic error
%skillf = @(x,y)((1 - (sum((x-y).^2)) / ... %index of agreement
%    (sum((x-abs(mean(y))).^2 + ...
%    (y-abs(mean(y))).^2))));
skillf = @(x,y)(1 -(sum((x-y).^2) / ...
      sum((abs(x-mean(y))+ abs(y-mean(y))).^2)));

%skillf = @(x,y)(1 -( sum((x-y).^2) / ...
%      sum(((abs(x-mean(y)))+ (abs(y-mean(y)))).^2)));
skill.d=skillf(pred,obs);

