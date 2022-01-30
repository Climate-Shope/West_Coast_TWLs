function [Out]=MorphLove_RegSlope(wave_profile,ESI,ESI_Array,Rey,OBJECTID,UTM);


%% Load in particulars
zzz=str2num(OBJECTID);
nn=cell2mat(ESI_Array(:,1));
mm=find(ismember(nn,zzz));
ESI_Better=mm(1);
%% Prefilter the profile before calculations

if length(wave_profile.depth(1,:))==1;
    wave_profile.depth=wave_profile.depth';
end

%Prefilter the profile to remove any inital decrease in slope, often caused
%by hydroprocessing of the DEM dataset
temp=diff(wave_profile.depth);
if temp(1)<0;
temp=find(temp>0);
temp=1:temp(1);
wave_profile.depth(temp)=[];
wave_profile.xpf(temp)=[];
wave_profile.ypf(temp)=[];
wave_profile.L=1:length(wave_profile.depth);
bass=1;
else
    temp=[];
end
if isempty(temp)
    addval=0;
else
    addval=length(temp);
end

MHW=wave_profile.mhw;
if wave_profile.mhw_Lpos>length(wave_profile.depth)
    MHW_loc=find(wave_profile.L==wave_profile.mhw_Lpos);%will change depending on profile
else
    MHW_loc=wave_profile.mhw_Lpos;
end
%If the previous section modified the profile, need to change the MHW
%location 
if exist('treb','var') || exist('bass','var')
    z=wave_profile.depth;
    x=1:length(z);
    h=repmat(MHW,1,length(z));
    
    L1=[x;z];
    L2=[x;h];
    P=InterX(L1,L2);
    if ~isempty(P);
    if P(1,1)>400;
        P=[];
    end
    end
    if ~isempty(P);
        P=round(P(1));
        MHW_loc=P;
    else
        MHW_loc=1;
    end
    
end

if MHW_loc<1;
    MHW_loc=1;
end
MSL=wave_profile.msl;

%check for riprap in front of the seawall
if ~isempty(Rey);
    if strcmp(ESI_Array(ESI_Better,5),'6B: Riprap');
        change_riprap=1;
    else
        change_riprap=0;
    end
end


%% Step 1, detrend and calculate the maxima (assuming this isnt a a dune
%environment), then calculate the detrend profile from 1 to 50 more than the
%calculated maxima. 

MHW=wave_profile.mhw;
depth1=wave_profile.depth;

[~,temp]=max(depth1);
[pks,locs] = findpeaks(wave_profile.depth,[],25,MHW);

depth1(temp:end)=depth1(temp);


%1. Find the most shoreward maximum

if isempty(Rey)
    %i.e. not a seawall focused profile
    %need to set up the profile smoothing to be dependent on  profile type
    %and unique threhsolds. These are determined by though testing by
    %region and will need to be modified by the user for different
    %locations
    if strcmp(ESI,'3A') || strcmp(ESI,'4') || strcmp(ESI,'5') && max(wave_profile.depth)<7%low lying sand to gravel beach
        result = DouglasPeucker([wave_profile.L;depth1],0.5);
    elseif strcmp(ESI,'3A') || strcmp(ESI,'4') || strcmp(ESI,'5')%Sand to gravel beach
        result = DouglasPeucker([wave_profile.L;depth1],1.0);
         elseif strcmp(ESI,'6B') || strcmp(ESI,'6A')%for engineered profiles
        result = DouglasPeucker([wave_profile.L;depth1],1.0); 
    elseif max(wave_profile.depth)>25 && length(wave_profile.depth)>100;
        result = DouglasPeucker([wave_profile.L;depth1],2.5);
    else
        
        result = DouglasPeucker([wave_profile.L;depth1],1.5);
    
    end
%%% find profile maximum 
if strcmp(ESI,'6B') && ~strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches');
    [pks,locs] = findpeaks(result(2,:),0.8,1,3.5);
elseif strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches');
    [pks,locs] = findpeaks(result(2,:),0.3,1,3.5);
elseif ~strcmp(ESI,'3A') && ~strcmp(ESI,'4')
    [pks,locs] = findpeaks(result(2,:),1,1,4); %Must have a 1 m prominence 
else
   [pks,locs] = findpeaks(result(2,:),0.8,1,3.5); %Must have a 1 m prominence  
end

% interpolate simplified profile to the same number of points as
% original profile
depth=interp1(result(1,:),result(2,:),wave_profile.L);
        
    %check to see if it is linear, if so, then no need to clip
    [p,S,mu] = polyfit(1:length(depth),depth,1);
    [y, delta] = polyval(p,1:length(depth),S,mu);
    te=sum(abs(depth-y));
    if te<1
        depth=wave_profile.depth;
        tiny=1;
    else
        tiny=0;
    end
    
    
    if ~isempty(locs)
    if strcmp(ESI,'3A') || strcmp(ESI,'4') || strcmp(ESI,'5')     
    locs=locs(1);
    locs=find(wave_profile.L==result(1,locs));
    else
    locs(result(2,locs)<10)=[];
    if ~isempty(locs)
        locs=locs(1);
        locs=find(wave_profile.L==result(1,locs));
    else
        locs=[];
    end  
    end
    end
    
 
    %If no locations of maxima are found and the profile is below 10m in
    %elevation, redoc the profile simplification and search for new
    %locations
    %
    if isempty(locs) && max(depth1)<10;
        if tiny==1;
        result = DouglasPeucker([wave_profile.L;depth1],0.01);
        [pks,locs] = findpeaks(result(2,:),[],1,MHW);
        else
            result = DouglasPeucker([wave_profile.L;depth1],0.5);
        
        [pks,locs] = findpeaks(result(2,:),1,1,MHW); %Must have a 1 m prominence 
        end
        depth=interp1(result(1,:),result(2,:),wave_profile.L);
 
         if ~isempty(locs)
            locs=locs(1);
            locs=find(wave_profile.L==result(1,locs));
         end
    end

    %% Add in bit to deal with very small maxima
    if strcmp(ESI,'3A') && ~isempty(locs) && max(depth1)>10;
        locs(wave_profile.depth(locs)<5)=[];
    end
    if strcmp(ESI,'3A') && ~isempty(locs);
        
        locs(wave_profile.depth(locs)<3.7)=[];
    end
    
    
if ~isempty(locs);
    maxi=locs(1);
else
    maxi=[];
end

%% Add in a bit to check for rocks in the way of the profile??
if ~isempty(maxi)
    
    if length(depth1)>50;
        yem=1;
    else
        yem=0;
    end
if yem == 1    
if depth1(maxi)<7 && depth1(50)>10
    maxi=[];
end
end
end
%%

%Next, limit the extrent of the profile such that we are always in a
%monotonically increasing scenario (ish)
if isempty(maxi)
    [tmp_max, tmp_ind]=max(depth);
    depth(tmp_ind:end)=tmp_max;

else
    depth(maxi:end)=wave_profile.depth(maxi);
end

%% Now calculate the cliff or dune top based on this edited profile

if isempty(maxi); 
    addx=(1:length(depth))+length(depth);
    %Create an artifical ramp to add to the end of the profile for use in
    %the crest identification
    ramp=wave_profile.depth(end)+(1:length(wave_profile.depth))*0.05;
    depth=[depth ramp];
else
    ramp=depth(end)+(1:length(depth))*0.05;
    depth=[depth ramp];
end

X=[1 length(depth)];
V=[depth(1) depth(end)];
Xq=[1:length(depth)];
Vq = interp1(X,V,Xq);
Values=depth-Vq;
[~,u]=max(Values);


if abs(wave_profile.depth(maxi)- wave_profile.depth(u))>1
    
    if u~=maxi && wave_profile.depth(u)<5 && mean(wave_profile.depth)>5
    maxi=maxi;
    elseif u~=maxi && wave_profile.depth(u)<5 && wave_profile.depth(u)<3.5 && mean(wave_profile.depth)>4.7 && strcmp(ESI,'3A')
    maxi=maxi;
    else
    maxi=u;
    end
end

%% If there is a seawall 
else
   Rey=Rey(1,:);
   tiny=0;
   if numel(num2str(round(Rey(:,3))))<3
      xx= wave_profile.xpf;
      yy= wave_profile.ypf;
      utmzone=repmat(wave_profile.UTMZONE,length(xx),1); 
       [Lat,Lon] = utm2deg(xx,yy,utmzone);
       xcomp=(Rey(:,2)-Lon).^2;
       ycomp=(Rey(:,3)-Lat).^2;
   else
       xcomp=(Rey(:,2)-wave_profile.xpf).^2;
       ycomp=(Rey(:,3)-wave_profile.ypf).^2;   
   end
   
   
    %Approximate the location along the profile of the seawall location 
    distance=sqrt(xcomp+ycomp);
    %define the maximum as the location of the seawall 
    [~,maxi]=min(abs(distance));
    maxi_org=maxi;
    
    %redo profile simplification and finding the maximum
    if maxi == 1 && maxi < MHW_loc
    result = DouglasPeucker([wave_profile.L;depth1],1.0);   
    [pks,locs] = findpeaks(result(2,:),0.5,1,3.5);   
     depth=interp1(result(1,:),result(2,:),wave_profile.L);
    if ~isempty(locs)
    locs=locs(1);
    locs=find(wave_profile.L==result(1,locs));
    end
    
     if isempty(locs) && max(depth1)<10;
        result = DouglasPeucker([wave_profile.L;depth1],0.5);
        [pks,locs] = findpeaks(result(2,:),1,1,MHW); %Must have a 1 m prominence
        depth=interp1(result(1,:),result(2,:),wave_profile.L);
         if ~isempty(locs)
            locs=locs(1);
            locs=find(wave_profile.L==result(1,locs));
         end
     end

    if ~isempty(locs);
        maxi=locs(1);
    else
        maxi=[];
    end  
           
    end
    
    if isempty(maxi)
        result = DouglasPeucker([wave_profile.L;depth1],0.3);
        [pks,locs] = findpeaks(result(2,:),0.5,1,MHW); %Must have a 1 m prominence 
        depth=interp1(result(1,:),result(2,:),wave_profile.L);
         if ~isempty(locs)
            locs=locs(1);
            locs=find(wave_profile.L==result(1,locs));
         end
          if ~isempty(locs);
            maxi=locs(1);
          else
            maxi=[];
          end  
         
    
    end
    
    %Do a local find peaks to grab the correct
    %maximum/crest by detrending the profile following Palaseanu-Lovejoy and other method  
    if ~isempty(maxi)
    %base detrend length on the location of the maximum location estimation
    if maxi+20<length(wave_profile.depth);
    tmp=detrend(wave_profile.depth(1:maxi+20));
    else
        tmp=detrend(wave_profile.depth(1:end));
    end
    tmp1=tmp;
    [pks,locs] = findpeaks(tmp,[],10,0.21);
    pks(locs<maxi-10)=[];
    pks(locs>maxi+20)=[];
    locs(locs<maxi-10)=[];
    locs(locs>maxi+20)=[];
    
    %try again if no peaks are found 
    if isempty(locs)
    [pks,locs] = findpeaks(tmp,[],5,0.21);
    pks(locs<maxi-7)=[];
    pks(locs>maxi+20)=[];
    locs(locs<maxi-7)=[];
    locs(locs>maxi+20)=[];   
    end

    if ~isempty(locs) 
        if ~isempty(maxi)
        %localize data to around the maxi estimation    
        pks(locs>maxi+20)=[];
        pks(locs<maxi-10)=[];
        locs(locs>maxi+20)=[];
        locs(locs<maxi-10)=[];
        xcomp=(locs-maxi).^2;
        ycomp=(locs-maxi).^2;
        distance=sqrt(xcomp+ycomp);
        if ~isempty(distance);
        [~,tmp]=min(abs(distance));
        maxi=locs(tmp);
        
        if tmp1(maxi)>2
            tmp2=tmp1;
            tmp2(1:maxi)=NaN;
            a=find(tmp2<2);
            if ~isempty(a)
            maxi=a(1);
            end
        elseif tmp1(maxi)>1
            maxi=maxi+1;
        end
        
        else
        maxi=[];
        
        Rey=[];
        
        end
        else
        maxi=locs(1);
        end
    

    else
        Rey2=maxi;
        maxi=[];
        Rey=[];
    end
    end
end
%% Now find maximum if maxi is empty

MHW=wave_profile.mhw;
depth1=wave_profile.depth;

%temp maximum location on profile
[~,temp]=max(depth1);
[pks,locs] = findpeaks(wave_profile.depth,[],25,MHW);

depth1(temp:end)=depth1(temp);


if isempty(maxi) || maxi ==1
    
    if exist('u','var')
        dess=wave_profile.depth(u);
        if dess>7
            dess=1;
        else
            dess=0;
        end
    else
        dess=0;
    end
    
    
    if mean(wave_profile.depth)<7 && dess~=1
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.2);
    else
    result = DouglasPeucker([wave_profile.L;depth1],1.0);    
    end
    
   %now find maximum potential maximum location based on profile type and parameters
   %as determined regionally by testing following much of the same
   %methodology    
if strcmp(ESI,'6B') && ~strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches');
    [pks,locs] = findpeaks(result(2,:),0.8,1,3.5);

elseif strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches') && mean(wave_profile.depth)<7 && dess~=1;
    [pks,locs] = findpeaks(result(2,:),0.2,1,3.25);
elseif strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches');
    [pks,locs] = findpeaks(result(2,:),0.3,1,3.5);
elseif mean(wave_profile.depth)<7
    [pks,locs] = findpeaks(result(2,:),0.3,1,3.5); %Must have a 1 m prominence 
else
    [pks,locs] = findpeaks(result(2,:),1,1,3.5); %Must have a 1 m prominence 
end
    depth=interp1(result(1,:),result(2,:),wave_profile.L);

 if ~isempty(locs)
    if strcmp(ESI,'3A') || strcmp(ESI,'4') || strcmp(ESI,'5')  
        
        if exist('Rey2','var')
            locf=[];
            for kkk=1:length(locs)
            locz=find(wave_profile.L==result(1,locs(kkk)));
            locf=[locf locz];
            end
            locs=locf;
            ymp=abs(locs-Rey2);
            [~,ymp]=min(ymp); locs=locs(ymp);
        else
        
            locs=locs(1);
            locs=find(wave_profile.L==result(1,locs));
        end
    
    else
    locs(result(2,locs)<10)=[];
    if ~isempty(locs)
        locs=locs(1);
        locs=find(wave_profile.L==result(1,locs));
    else
        locs=[];
    end  
    end
 end
    
 

    if isempty(locs) && max(depth1)<10;
        if exist('dess','var');
        if mean(wave_profile.depth)<7 && dess~=1
            result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.2);
            [pks,locs] = findpeaks(result(2,:),0.5,1,MHW);
        else
        result = DouglasPeucker([wave_profile.L;depth1],0.5);
        [pks,locs] = findpeaks(result(2,:),1,1,MHW);
        end
        else
        [pks,locs] = findpeaks(result(2,:),1,1,MHW); %Must have a 1 m prominence 
        end
        depth=interp1(result(1,:),result(2,:),wave_profile.L);
         if ~isempty(locs)
             
            locs=locs(1);
            locs=find(wave_profile.L==result(1,locs));
         end
    end
    
     %% Add in bit to deal with very small maxima
    if strcmp(ESI,'3A') && ~isempty(locs) && max(depth1)>10 && mean(wave_profile.depth)>7;
        locs(wave_profile.depth(locs)<5)=[];
    end
     if strcmp(ESI,'3A') && ~isempty(locs);
        
        locs(wave_profile.depth(locs)<3.7)=[];
    end

if ~isempty(locs);
    maxi=locs(1);
else
    maxi=[];
end

%% Check for rocks in the way of the profile
if ~isempty(maxi)
    if length(depth1)>50;
        yem=1;
    else
        yem=0;
    end
if yem == 1    
if depth1(maxi)<7 && depth1(50)>10
    maxi=[];
end
end
end

%Next, limit the extrent of the porfile such that we are always in a
%monotonically increasing (ish)
if isempty(maxi)
    [tmp_max, tmp_ind]=max(depth);
    depth(tmp_ind:end)=tmp_max;

else
    depth(maxi:end)=wave_profile.depth(maxi);
end

%Now calculate the cliff or dune top based on this edited profile


if isempty(maxi); 
    if max(wave_profile.depth)>3.25
     addx=(1:length(depth))+length(depth);
     
     if exist('u','var')
     if strcmp(ESI,'3A') && u<20 && wave_profile.depth(u)<3.2
         ramp=wave_profile.depth(end)+(1:length(wave_profile.depth))*0.02;
     else
         ramp=wave_profile.depth(end)+(1:length(wave_profile.depth))*0.03;
     end
     else
        ramp=wave_profile.depth(end)+(1:length(wave_profile.depth))*0.03; 
     end
     
     if exist('dess','var');
        if mean(wave_profile.depth)<7 && dess~=1 && length(wave_profile.depth)>300
            ramp=[];
        end
     end

    depth=[depth ramp];
    end
end
    

% %maxi=maxi+1;

if ~isempty(Rey) && isempty(maxi) 
    X=[1 length(depth)];
    V=[depth(1) depth(end)];
    Xq=[1:length(depth)];
    Vq = interp1(X,V,Xq);
    Values=depth-Vq;
    [pks,locs] = findpeaks(Values,[],1,-10);
    
    if ~isempty(locs)
        [~,rrr]=min(abs((maxi_org-locs)));
        maxi=locs(rrr);
    else
        maxi=maxi_org;
    end

elseif   isempty(Rey) && exist('maxi_org','var')
    [pks,locs] = findpeaks(depth,[],1,3);
    
    if ~isempty(locs)
        ddd=abs(locs-maxi_org);
        locs(ddd>20)=[];   
    end
    
    
    if ~isempty(locs);
        [~,ddd]=min(abs(locs-maxi_org));
        maxi=locs(ddd);
    else
       X=[1 length(depth)];
        V=[depth(1) depth(end)];
        Xq=[1:length(depth)];
        Vq = interp1(X,V,Xq);
        Values=depth-Vq;
        
        if exist('Rey2','var')
           [pks,locs] = findpeaks(Values,[],10,0.5);
           if ~isempty(locs) && ~isempty(maxi)
           glo=abs(locs-maxi);
           [~,tin]=min(glo);
           maxi=locs(tin);
           else
            [~,maxi]=max(Values);    
           end
        else
             [~,maxi]=max(Values); 
        end
        
       
    end
    
    
elseif max(wave_profile.depth)>3.25
X=[1 length(depth)];
V=[depth(1) depth(end)];
Xq=[1:length(depth)];
Vq = interp1(X,V,Xq);
Values=depth-Vq;

if strcmp(ESI,'6B')
cut_idx=find(depth<3.25);

if length(cut_idx)~=length(depth)
    Values(cut_idx)=NaN;
end
    


elseif strcmp(ESI,'3A') && max(wave_profile.depth)>3.7
    cut_idx=find(depth<3.7);Values(cut_idx)=NaN;
end


if strcmp(ESI,'6B')
    [pks,locs] = findpeaks(Values,[],1,1);
    if ~isempty(locs)
        locs(wave_profile.depth(locs)<3.25)=[];
        if ~isempty(locs)
           maxi=locs(1);
        else
            [~,maxi]=max(Values);
        end
    else
        [~,maxi]=max(Values);
    end
else
    [~,maxi]=max(Values);
end


if maxi==1 && exist('u','var');
    maxi=u;
end
    
else
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pks,locs] = findpeaks(wave_profile.depth,[],1,1);
    
    if~isempty(locs)
    [~,yy]=max(pks);
    maxi=locs(yy);
    else
        maxi=1;
    end

end
end

if maxi>length(wave_profile.depth) && exist('u','var')
    maxi=u; 
end


%% Look for a more seaward maximum if Rey (seawall present) exists
if ~isempty(Rey);
[pks,locs] = findpeaks(wave_profile.depth,[],5,1);
if ~isempty(locs)
    pks(locs>maxi)=[];
    locs(locs>maxi)=[];
    
    locs(pks<wave_profile.depth(maxi)+0.8)=[];
    pks(pks<wave_profile.depth(maxi)+0.8)=[];
    
    if~isempty(locs)
        [~,zz]=max(pks);
        maxi=locs(zz);
        
    end

end

end

%% Now to select the correct toe location
%limit the profile length to 2x the distance to maxi, but then I want to test a
%few and find a convergence point
depth2=wave_profile.depth;

if maxi>70 && wave_profile.depth(maxi)>20
ramp=wave_profile.depth(maxi)+([1:length(maxi+1:length(wave_profile.depth))])*0.04;
else
ramp=wave_profile.depth(maxi)+([1:length(maxi+1:length(wave_profile.depth))])*0.02;
end
depth2(maxi+1:end)=ramp;

%simplify the profile a little 1st
if strcmp(ESI,'3A') || strcmp(ESI,'5')
result = DouglasPeucker([wave_profile.L;depth2],0.2);
elseif ~isempty(Rey)
result = DouglasPeucker([wave_profile.L;depth2],0.1);
else
result = DouglasPeucker([wave_profile.L;depth2],0.6);    
end

depth2=interp1(result(1,:),result(2,:),wave_profile.L);

if maxi<5 && wave_profile.depth(1)>3.5
    toe=1;
    toez=1;
elseif maxi==1 && strcmp(ESI,'3A') || tiny==1;
    toe=1;
    toez=1;
else   
    %define endpoints for use in an iterative approach to define the toe
    %location 
    endpos=[maxi*2;round(maxi*1.8);round(maxi*1.6);round(maxi*1.4);round(maxi*1.2);maxi];

    toe=[];
for ii=1:length(endpos);

gum=endpos(ii); 
if gum > length(wave_profile.depth)
   gum=length(wave_profile.depth);
end
newdepth=depth2(1:gum);
%Find all locations on profile that are above the defined thresholds
%(determined through testing)

if strcmp(ESI,'5')
    flump=find(newdepth>10);
elseif strcmp(ESI,'6A') || strcmp(ESI,'6D')
    flump=find(newdepth>12);
elseif strcmp(ESI,'6B') 
    flump=find(newdepth>5.5);
    if ~isempty(flump)
    fixslope=movingslope(newdepth,3);
    fixslope=movingslope(fixslope,3);
    [~,cat]= findpeaks(fixslope,[],5,0);
    
    if ~isempty(cat)
        if cat(1)>flump(1)
            if (newdepth(cat(1))+1)>max(newdepth);
                flump=[];
            else
               flump=find(newdepth>(newdepth(cat(1))+1));  
            end
           
        end
        
        if ~isempty(flump)
        if newdepth(flump(1))>10
            flump=find(newdepth>7.5); 
        end
        end
        
    end
    end
else
flump=find(newdepth>8);
end


X=[1 length(newdepth)];
V=[newdepth(1) newdepth(end)];
Xq=[1:length(newdepth)];
Vq = interp1(X,V,Xq);
Values=newdepth-Vq;
Values(maxi+1:end)=[];


if ~isempty(flump)
    if flump(1)<maxi+1
       Values(flump:end)=[]; 
    end    
end
    
Valin=Values*-1;
[~,Valin]= findpeaks(Valin,[],10,-4);

%Fix so the toe does not exceed 10m 

if strcmp(ESI,'5')
    vmp=find(newdepth>13);
else
vmp=find(newdepth>10);
end
if ~isempty(vmp)
    cmp=vmp(1);
    if cmp<length(Values);
     Values=Values(1:vmp(1));
    end
end


if strcmp(ESI,'3A') || strcmp(ESI,'4');

    Values=Values;
    a=find(wave_profile.depth<2.5);
    
    %% For Ventura beaches in particular %%
    
    if ~isempty(a);
    ab=diff(a);
    b=find(ab>1);
    
    if isempty(b);
        b=a(end);
    end
    
    if ~isempty(b)
    b=b(1)-1;
    Values(1:b)=NaN;
    end
        a=find(wave_profile.depth>10);
        if ~isempty(a);
            a=a(1);
            Values(a:end)=NaN;
        end

    
    a=find(wave_profile.depth>wave_profile.depth(maxi)-1);
    if ~isempty(a);
        a=a(1);
        Values(a:end)=NaN;
    end
    
    end
end
    
if maxi==length(depth1);
    a=find(Values>1.5);
    if ~isempty(a);
         a=a(1);
        Values(a:end)=NaN;
    end
end

[ump,guess]=min(Values);


    Valin(Valin>guess)=[];
if~isempty(Valin)    
if ~ismember(guess,Valin)
    guess=Valin(end); %Pick closest inflection
end

elseif isempty(Valin) && ~isempty(Rey)
    ff=movingslope(newdepth,3);ff=movingslope(ff,3);
    [~,mess]= findpeaks(ff,[],10,0);
    
    if ~isempty(mess)
        mess(mess>maxi)=[];
        
        if ~isempty(mess)
            guess=mess(end);
        end
        
    end
    
end


toe=[toe;guess];

end


%Determine best toe to use of possible locations
toe(toe==1)=[];
if strcmp(ESI,'1B') && maxi>175 && ~strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches');
    toe(toe>125)=[];
end

toez=toe;

if isempty(toe);
    der=diff(depth2);
    der2=diff(depth2);
end

if~isempty(toe);  
    if strcmp(ESI,'3A') && length(toe)>1 || strcmp(ESI,'4') && length(toe)>1
    toe(1)=[];
    end
end

if~isempty(toe); 
  ggg=movingslope(depth2,3); hhh=movingslope(ggg,3);
  if depth2(maxi)<8 && strcmp(ESI,'3A') || depth2(maxi)<8 && strcmp(ESI,'4')
      [~,locs] = findpeaks(hhh,[],1,0.025);
  else
  [~,locs] = findpeaks(hhh,[],1,0.04);
  end
  if~isempty(locs); 
    toe=toe(ismember(toe,locs));
    
    if isempty(toe) && length(toez)>1
    toe=toez(2:end);
    end
  end
end

if~isempty(toe);  
    
    if maxi > 50 && wave_profile.depth(maxi)>15 && strcmp(ESI,'3A') 
        toe=toe(end);
    else
        toe=toe(1);
    end

else
    toe=1;
end
 
end
if toe<MHW_loc && MHW < 10;
    toe=MHW_loc;
end

 %Now if there is a seawall, and there isnt riprap, 
 %change the toe location to something that is closer to the seawall
 %location, with a max of 3m in front of the seawall crest location
 
 if ~isempty(Rey);
     if change_riprap==0; %Should not need to modify if riprap in front, code should pick it up
         depth3=depth2;
         depth3(maxi+1:end)=[];

         depth3(1:maxi-20)=NaN; %Limit the search radius around the maximum
         a=diff(depth3);%1st derivative
         gg=atan(a)*(180/pi);
         a=diff(a);%2nd derivative 
         [~,locs] = findpeaks(a,[],3,0.06);
         
         
       if max(gg)>3;
         
         if ~isempty(locs);
         locs=locs+1;%correct for the nature of diff function
         end
         if numel(locs)==1;
             toe=locs;
         elseif numel(locs)>1;
             toe=locs(end); %Selects the one closer to the crest of seawall
         elseif isempty(locs);
             toe=maxi-3;
         end
       end
     end
 end
         
         

%for riprap
if strcmp(ESI,'6B') && ~strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches');
result = DouglasPeucker([wave_profile.L;depth1],0.1);
fepth=interp1(result(1,:),result(2,:),wave_profile.L); 
a=diff(fepth);
a=diff(a);
[pks,locs] = findpeaks(a,[],1,0.2);

locs(locs>toe)=[];
if ~isempty(locs);
toe=locs(1);
end

%now with the new toe, find a new maximum for the riprap
result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.1);
fepth=interp1(result(1,:),result(2,:),wave_profile.L); 
[pks,locs] = findpeaks(fepth,[],1,3.5);
locs(locs<toe)=[];
if ~isempty(locs) && wave_profile.depth(maxi)>15
maxi=locs(1);
end



end

%% Now add in stockdon toe and similar maximum
%the toe here is defined as the first maximum in the second derivative
%meeting the conditions delow, as determined by testing for the region
result = DouglasPeucker([wave_profile.L;depth1],0.2);
fepth=interp1(result(1,:),result(2,:),wave_profile.L); 
a=diff(fepth);
b=diff(a);%2nd derivative approximation

if strcmp(ESI,'3A') && toe>100 || strcmp(ESI,'4') && toe>100
[pks,locs] = findpeaks(b,[],1,0.07);
locs(locs<30)=[];
else
[pks,locs] = findpeaks(b,[],1,0.15);    
end

a=find(a<0);
locs(ismember(locs,a))=[];


locs(locs>toe)=[];
locs(locs<MHW_loc)=[];
locs(locs<2.5)=[];
if ~isempty(locs);
toe_onshore=locs(1);
else
   toe_onshore=toe; 
end

%now with the new toe, lets find a new maximum 
result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.1);
fepth=interp1(result(1,:),result(2,:),wave_profile.L); 
[pks,locs] = findpeaks(fepth,[],1,3.5);
if exist('toe_onshore','var')
if~isempty(toe_onshore);
locs(locs<toe_onshore)=[];
else
  locs(locs<toe)=[];  
end
end
if ~isempty(locs)
maxi_onshore=locs(1);
end

%% Now there are two toe estimates, one from the modified Palasenu-Lovejoy 
%method and one fform the first significant onshore inflection point
%this section accounts for the toe onshore and toe discrepancy 
if toe-toe_onshore > 60 && toe_onshore<toe && strcmp(ESI,'1B');
    AT=atan(movingslope(depth1,3))*(180/pi);
    AT=AT(toe_onshore:toe);
    if max(AT>17)
      toe=toe_onshore;
    end
end

%2A profiles provide a more complex morphology for the toe
if strcmp(ESI,'2A') && toe>10;
    
    %First check toe location
    depth5=wave_profile.depth;
    result = DouglasPeucker([wave_profile.L;depth5],0.1);
    depth5=interp1(result(1,:),result(2,:),wave_profile.L);
    depth5=depth5(1:toe);
    
    
    
    a=diff(depth5);
    b=diff(a);
    
    if MHW_loc<=length(depth5)
     c=find(depth5<depth5(MHW_loc));
  
    b=find(b>0.3);
    b=b+1;
     b(ismember(b,c))=[];
    else
        b=[];
    end
   
    
    if ~isempty(b);
        
        toe=b(end);
    end
    
    
    depth5=wave_profile.depth;
    result = DouglasPeucker([wave_profile.L;depth5],0.5);
    depth5=interp1(result(1,:),result(2,:),wave_profile.L);
    depth5=depth5(1:toe);
    tmp=detrend(depth5);
    [~,tmp]=max(tmp);
      bench_loc=tmp;
    bermwidth=toe-bench_loc;
    bench_toe_loc=wave_profile.mhw_Lpos;
    if bench_toe_loc<0;
        bench_toe_loc=1;
    end
    
    if bench_loc==toe;
        bench_loc=bench_toe_loc;
         bermwidth=toe-bench_loc;
    end
    
    
    %%%Now add a check to ensure that the berm is roughly less than 1/15
    %%%slope as somewhat specified by the TAW guidelines
    if bench_loc>0;
        Clope=(wave_profile.depth(toe)-wave_profile.depth(bench_loc))/(bermwidth);
    else
        Clope=(wave_profile.depth(toe)-wave_profile.depth(1))/(toe-1);
    end
    
    %If the slope is not reasonaly horizontal, then do not consider it a berm condition roughly from the TAW guidelines
    if abs(Clope)>(1/12.5);
        bench_loc=[]; 
        bermwidth=[];
     end
    
    if toe==1; %There cant really be a berm to measure in this case
        bench_loc=[]; 
        bermwidth=[]; 
        bench_toe_loc=[];
    end
    
    if ~isempty(bench_toe_loc) && isempty(bench_loc) && bench_toe_loc > 0;
        bench_loc=bench_toe_loc;
        bermwidth=toe-bench_loc; 
        Clope=(wave_profile.depth(toe)-wave_profile.depth(bench_loc))/(bermwidth);
    if Clope>(1/15) || bermwidth<2; %If the slope is not reasonaly horizontal, then do not consider it a berm condition roughly from the TAW guidelines
        bench_loc=[]; 
        bermwidth=[];
    end
    end
    
    if ~isempty(toe)
    if bench_loc>=toe
        bench_loc=[]; 
        bermwidth=[];
        bench_toe_loc=[];
    end
    end
    
else
    bench_loc=[];
    bermwidth=[];
    bench_toe_loc=[];
end




%One more check for Rocky Coastlines 
if strcmp(ESI,'1A') || strcmp(ESI,'1B')|| strcmp(ESI,'2A') || strcmp(ESI,'6D') || strcmp(ESI,'6B') || strcmp(ESI,'6A')
    
    if toe<10 && toe > 1
    aa=wave_profile.depth(1:toe-1);
    aa=find((diff(aa))>0.67);%determined through testing
    if ~isempty(aa);
        toe=1;
    end
    end
end


%% Replace toe location with toe onshore if there is a length where buildings have been removed from the DEM
% or the calculated toe location is very far onshore and the estimates
% differ significantly


if toe>100 && toe-toe_onshore>50 && strcmp(ESI,'6B')
result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.5);
jepth=interp1(result(1,:),result(2,:),wave_profile.L); 
jlope=movingslope(jepth,3);
jlope(toe:end)=[];

if abs(mean(jlope(end-50:end)))<0.01
    toe=toe_onshore;
end

elseif strcmp(ESI,'6B') && toe-toe_onshore>30 || strcmp(ESI,'1B') && toe-toe_onshore>30
    if strcmp(ESI_Array{str2num(OBJECTID),6},'revetment')
        toe=toe_onshore;
    end
 
elseif strcmp(ESI,'3A') && toe>90 && toe>toe_onshore|| strcmp(ESI,'4') && toe>90 && toe>toe_onshore
    toe=toe_onshore;
    mod_toe=1;   
    
elseif toe-toe_onshore>35 && toe_onshore>40 && strcmp(ESI,'4')
    toe=toe_onshore;
    
elseif strcmp(ESI,'5') && toe-toe_onshore>20 && wave_profile.depth(toe)-wave_profile.depth(toe_onshore)>3 && max(wave_profile.depth(toe_onshore:toe))>wave_profile.depth(toe)
    toe=toe_onshore;
    

end

%% Next up, utilize the ESI and the slopes around the toe to determine
% the slopes that will be needed for runup method determination 


%% Look for the cliff junction (the cliff crest) if there is one
%all threhsolds have been determined for the region through testing
if maxi<5 && wave_profile.depth(1)>3.5
     cliffjun_loc=maxi;
elseif maxi==1 && strcmp(ESI,'3A')
    cliffjun_loc=maxi;
elseif length(wave_profile.depth)<10 && maxi<5 && wave_profile.depth(maxi)<3.5
    cliffjun_loc=maxi;
elseif maxi<4
    cliffjun_loc=maxi;    
else
if maxi == length(wave_profile.depth);
depth4=wave_profile.depth(1:maxi);
else
 depth4=wave_profile.depth(1:maxi+1);   
end

X=[1 length(depth4)];
V=[depth4(1) depth4(end)];
Xq=[1:length(depth4)];
Vq = interp1(X,V,Xq);
Values=depth4-Vq;
Values(1:toe-1)=NaN;
Values(end-3:end)=NaN;

if strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches') && toe>100;
    [~,i] = findpeaks(Values,[],5,-20);
    if ~isempty(i)
    i=i(1);
    else
     [~,i]=max(Values);
    end
elseif strcmp(ESI,'1B') && strcmp(ESI_Array(ESI_Better,5),'5: Mixed sand and gravel beaches') && toe>30;
    [~,i] = findpeaks(Values,[],5,-20);
    if ~isempty(i)
    i=i(1);
    else
     [~,i]=max(Values);
    end    
else
    [~,i]=max(Values);
end

if i==maxi && length(depth4)-i<5
    
    depth10=wave_profile.depth(toe:maxi);  
    X=[1 length(depth10)];
    V=[depth10(1) depth10(end)];
    Xq=[1:length(depth10)];
    Vq = interp1(X,V,Xq);
    Values=depth10-Vq;
    [~,i]=max(Values);
    i=i+toe;
if maxi<5 && wave_profile.depth(1)>3.5
     cliffjun_loc=maxi;
else
if maxi == length(wave_profile.depth);
depth4=wave_profile.depth(1:maxi);
else
 depth4=wave_profile.depth(1:maxi+1);   
end

X=[1 length(depth4)];
V=[depth4(1) depth4(end)];
Xq=[1:length(depth4)];
Vq = interp1(X,V,Xq);
Values=depth4-Vq;
Values(1:toe-1)=NaN;
    
end
end

bbb=find(Values>0);
%Calculate general slope between maximum and toe location
if toe<1
dddd=(wave_profile.depth(maxi)-wave_profile.depth(1))/(maxi-1);
else
dddd=(wave_profile.depth(maxi)-wave_profile.depth(toe))/(maxi-toe);   
end
dddd=atan(dddd)*(180/pi);
if isempty(bbb) && strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches') && toe>100;
    [pks,locs] = findpeaks(Values,[],20,-20);
elseif toe>50 && strcmp(ESI,'4') && maxi-toe>50 && dddd<20;
    [pks,locs] = findpeaks(Values,[],20,-15);
    zipzip=1;
elseif isempty(bbb) && strcmp(ESI,'4')
    [pks,locs] = findpeaks(Values,[],20,-4);
elseif isempty(bbb);
    [pks,locs] = findpeaks(Values,[],20,-4);
else
[pks,locs] = findpeaks(Values,[],20,-3);
end

locs(locs>i)=[];
locs(depth4(locs)<4)=[];
if ~isempty(locs)
    if length(locs)>1 
        i=locs(end-1);
    elseif exist('zipzip','var')
        i=locs(1);
    end
end

val=Values(i);
cut=val-2;
if val<0
[pks,locs] = findpeaks(Values,[],10,cut);
else
    if cut<0
        cut=0;
    end
    [pks,locs] = findpeaks(Values,[],10,cut);
end

if~isempty(locs)
    pks(locs>maxi)=[];
    pks((val-Values(locs))>1.5)=[];
    locs(locs>maxi)=[];
    locs((val-Values(locs))>1.5)=[];
    
      
    if~isempty(locs)
        if max(pks)-pks(1)<7 && max(pks)-pks(1)~=0
           [red,green] = findpeaks(Values,0.1,10,-5); 
            blue=find(ismember(locs,green));
            blue=locs(blue);
            
            if ~isempty(blue);
                i=blue(1);
            else
                i=locs(1);
            end
        elseif strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches') && toe>100
            i=locs(1);
        else
            [~,pppp]=max(pks);
            i=locs(pppp);
        end
    end
end

if strcmp(ESI,'6A');
[red,green] = findpeaks(Values,0.1,10,3);    
elseif strcmp(ESI,'3A') && wave_profile.depth(maxi)>25 && maxi-toe<40;
    [red,green] = findpeaks(Values,[],2,-6);
else
    [red,green] = findpeaks(Values,0.1,10,-6);
end


if strcmp(ESI,'2A') && toe<10 && isempty(bench_loc);
    green(wave_profile.depth(green)<5)=[];
end
if ~isempty(green)
    
    cliffjun_loc=green(1);
else
cliffjun_loc=i;
end


if wave_profile.depth(cliffjun_loc)<0.15
    cliffjun_loc=maxi;
end

end


if maxi-toe<5 && ~isempty(Rey) && cliffjun_loc<toe
    cliffjun_loc=maxi; 
end
    
    
%% Now look at the scenario where there may be another plausible toe
%  between the first pass toe and the cliff junction locations
%  toe_up is produced to be reviewed in script while debugging/evaluating 

result = DouglasPeucker([wave_profile.L;depth1],0.1);
fepth=interp1(result(1,:),result(2,:),wave_profile.L); 
a=diff(fepth);
a=diff(a);
[pks,locs] = findpeaks(a,[],1,0.15);

locs(locs<toe)=[];
locs(locs<MHW_loc)=[];
locs(locs>cliffjun_loc)=[];
if ~isempty(locs);
toe_up=locs(1);
else
   toe_up=toe; 
end

%% Create a catch for the toe in a number of places where there is the 
%potential for a platform in front of the estimated toe location
%As always, simplify profile first
if ~isempty(toe)
    depth6=wave_profile.depth;
    result = DouglasPeucker([wave_profile.L;depth6],0.3);
    depth6=interp1(result(1,:),result(2,:),wave_profile.L);
    depth6=depth6(1:toe);
    a=diff(depth6);
    c=find(a>0.4);
    if ~isempty(c);
    b=diff(a);
    [pks,locs] = findpeaks(b,[],[],0.25);
    locs=locs+1;
    locs(locs<MHW_loc)=[];
    
    if ~isempty(locs)
        extra_toe=locs(end);
    else
       extra_toe=[]; 
    end
    else
     extra_toe=[];
    end
else
    extra_toe=[];
end


%% Fix the MHW location in regards to new toe locations 
if MHW_loc>10 && toe>MHW_loc;
    L1=[1:length(wave_profile.depth);wave_profile.depth];
    L2=[[1 length(wave_profile.depth)];[wave_profile.mhw wave_profile.mhw]];
    P = InterX(L1,L2); 
    if ~isempty(P);
        P=round(P(1,:));
        P(P>toe)=[];
        MHW_loc=P(end);
    end       
end


%% Now to correct for a large slope before the toe
if ~strcmp(ESI,'2A') && isempty(Rey);%again Rey is a flag for seawalls
    if ~isempty(toe)
        if toe>5 && toe>MHW_loc
            a=wave_profile.depth(1:toe-1);
            b=movingslope(a,3);
            c=find(b>0.67);
            if~isempty(c)
                if c(1)<10
                   toe=MHW_loc;
                   toe_change=1;
                elseif c(1)< toe-3;
                   a=wave_profile.depth(1:c(1)+3);
                   b=diff(a); 
                   e=diff(b);
                   [pks,locs] = findpeaks(e,[],[],0.25);
                   locs(locs>c(1)+1)=[];
                   if~isempty(locs);
                       toe=locs(end);
                       toe_change=1;
                   else
                       toe=MHW_loc;
                   end
                
                   %If the toe changes, need to looks at the cliff junction
                   if exist('toe_change','var')
                     X=[1 length(depth4)];
                    V=[depth4(1) depth4(end)];
                     Xq=[1:length(depth4)];
                    Vq = interp1(X,V,Xq);
                    Values=depth4-Vq;
                    Values(1:toe-1)=NaN;
                    Values(end-3:end)=NaN;
                    [pks,locs] = findpeaks(Values,[],10,-4);
                    
                    if~isempty(locs)
                        [~,ccc]=max(pks);
                        cliffjun_loc=locs(ccc);
                    end
                    
                   end

                end 
                
                if exist('toe_change','var')
                     X=[1 length(depth4)];
                    V=[depth4(1) depth4(end)];
                     Xq=[1:length(depth4)];
                    Vq = interp1(X,V,Xq);
                    Values=depth4-Vq;
                    Values(1:toe-1)=NaN;
                    Values(end-3:end)=NaN;
                    [pks,locs] = findpeaks(Values,[],10,-4);
                    
                    if~isempty(locs)
                        [~,ccc]=max(pks);
                        cliffjun_loc=locs(ccc);
                    end
                    
                 end
                
                
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fix maxima that could be a rock at sea and not a crest
%% Find infleciton points

result = DouglasPeucker([wave_profile.L;depth1],1.0);
depth=interp1(result(1,:),result(2,:),wave_profile.L);

a=diff(depth);
b=diff(a);

[pks,locs] = findpeaks(b,0.01,1,0.1);
locs(locs<maxi-1)=[];
if ~isempty(locs)
inflect1=locs(1)+1;

inflect1_slope=(depth(inflect1-1)-depth(inflect1))/(((inflect1-1)-(inflect1)));

c=b*-1;
[pks,locs] = findpeaks(c,0.01,1,0.1);
locs(locs<inflect1-1)=[];
if ~isempty(locs) && inflect1_slope<0;
inflect2=locs(1)+1;
inflect2_slope=(depth(inflect2-1)-depth(inflect2+1))/(((inflect2-1)-(inflect2+1)));

inflect_dist=inflect2-inflect1;

%% 
if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4') && ~strcmp(ESI,'6B') && isempty(Rey)%no big rocks in the way in riprap the rocks are pre-intsalled
if inflect_dist<15 && min(wave_profile.depth(inflect1:inflect2))<7 && inflect2_slope>0 && max(depth1)>7 && inflect1<50; %Basically assume that the rock in front is part of the overall morphology
    
    %now search for maximum using the old methodology 
    ramp=wave_profile.depth(end)+(1:length(depth))*0.05;
    depth=[depth ramp];
    X=[1 length(depth)];
    V=[depth(1) depth(end)];
    Xq=[1:length(depth)];
    Vq = interp1(X,V,Xq);
    Values=depth-Vq;
    Values(1:inflect2)=NaN;
    Values(length(wave_profile.depth)+1:end)=NaN;
    [~,maxi]=max(Values);
    cliffjun_loc=maxi;
    twitch=1;
    
end
end
end
end
%Redo the cliff junction analysis with a new maximum location
if maxi==length(wave_profile.depth);   

depth4=wave_profile.depth(1:maxi);

if toe>toe_onshore || wave_profile.depth(toe_onshore)>8;
    flop=toe;
else
    flop=toe_onshore;
end

X=[flop length(depth4)];
V=[depth4(flop) depth4(end)];
Xq=[flop:length(depth4)];
Vq = interp1(X,V,Xq);
Values=depth4(flop:end)-Vq;

if strcmp(ESI,'2A')
    fff=150-flop;
    if fff>1
    Values(fff:end)=[];
    end
    tw=1;
end

if exist('inflect2','var');
    Values(1:inflect2+1)=NaN;
end

[~,i]=max(Values);
val=Values(i);
if exist('tw','var');
[pks,locs] = findpeaks(Values,[],10,1.5);

locs(wave_profile.depth(locs+flop)<5)=[];

else
 [pks,locs] = findpeaks(Values,[],10,5);   
end

if~isempty(locs)
    locs(locs>maxi)=[];
    locs((val-Values(locs))>1)=[];
    if~isempty(locs)
        i=locs(1);
    end
end


if Values(i)<1
    result = DouglasPeucker([wave_profile.L;depth4],1.0);
    depth7=interp1(result(1,:),result(2,:),wave_profile.L);
    a=diff(depth7);
    a=diff(a);
    a=a*-1;
    if strcmp(ESI,'3A')
        [pks,locs] = findpeaks(a,0.01,1,0.35);
    else
    [pks,locs] = findpeaks(a,0.01,1,0.5);
    locs((locs-toe)>70)=[];
    
    if isempty(locs)
      [pks,locs] = findpeaks(a,0.01,1,0.15);  
    end
    
    end
    locs=locs+1;
    locs(locs<flop)=[];
    if exist('inflect2','var');
    locs(locs<inflect2+1)=[];
    end
    
    
    if ~isempty(locs)
    locs=locs(1)+1; %Want the point just shoreward of the junction
    cliffjun_loc=locs;
    elseif exist('blue','var');
        if ~isempty(blue)
        cliffjun_loc=cliffjun_loc;
        else
          cliffjun_loc=maxi;  
        end
        
    elseif strcmp(ESI,'6B') && strcmp(ESI_Array(ESI_Better,5),'3A: Fine- to medium-grained sand beaches') && toe>100
        cliffjun_loc=cliffjun_loc;
    elseif exist('zipzip','var')
        cliffjun_loc=cliffjun_loc;
    else
     cliffjun_loc=maxi;
    end
else

if (i+flop)-cliffjun_loc<20 && (i+flop)-cliffjun_loc>-3
    cliffjun_loc=i+flop;

elseif strcmp(ESI,'3A') && (i+flop)-cliffjun_loc<20 && (i+flop)-cliffjun_loc>-15 && wave_profile.depth(i+flop)>9 
    zlope=atan((wave_profile.depth(i+flop)-wave_profile.depth(flop))/(i+flop-flop))*(180/pi);
    
    if zlope>45
       cliffjun_loc=i+flop;        
    end
    
end

if wave_profile.depth(cliffjun_loc)<0.15
    cliffjun_loc=maxi;
end


end 

result = DouglasPeucker([wave_profile.L;depth4],1.0);
depth7=interp1(result(1,:),result(2,:),wave_profile.L);
Dvec = movingslope(depth7,5);ff1=atan(Dvec)*(180/pi);
gg=diff(depth7);ff=atan(gg)*(180/pi);
hh=diff(ff);

if cliffjun_loc~=maxi && maxi-cliffjun_loc>3

uu=cliffjun_loc-1;
vv=hh*-1;
[pks,locs] = findpeaks(vv,0.01,1,5);
locs(locs>uu+2)=[];
locs(depth7(locs)<5)=[];
    if exist('inflect2','var');
    locs(locs<inflect2)=[];
    end

end
end

if wave_profile.depth(cliffjun_loc)<9 && ~strcmp(ESI,'3A') && ~strcmp(ESI,'4') && length(wave_profile.depth)<100 && exist('depth4','var')
    
    X=[1 length(depth4)];
    V=[depth4(1) depth4(end)];
    Xq=[1:length(depth4)];
    Vq = interp1(X,V,Xq);
    Values=depth4-Vq;
    Values(1:toe-1)=NaN;
    Values(end-3:end)=NaN;
    [~,i]=max(Values);
    
    i(i>maxi)=[];
    i(i<toe)=[];
    
    if ~isempty(i) && strcmp(ESI,'1A') || ~isempty(i) && strcmp(ESI,'2A');
        cliffjun_loc=i;
    else
        cliffjun_loc=maxi;
    end
        
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Addendum rock removal section 4/25/2019
     if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4')
         avg=2;
       L1=[[1:length(wave_profile.depth)];wave_profile.depth];
       L2=[[1 length(wave_profile.depth)];[avg avg]];
        P = InterX(L1,L2);
        
        if ~isempty(P);
            T=round(P(1,:)); 
            T(T>cliffjun_loc)=[];
        end
        
        if exist('T','var')
            if ~isempty(T)
                if length(T)>1
                T=T(end);
                if toe<T
                  if ~isempty(toez);
                     toez(toez<T)=[]; toez(toez>cliffjun_loc)=[];
                     if ~isempty(toez);
                     toe=toez(1);
                     toe_onshore=toe;
                     else
                     result = DouglasPeucker([wave_profile.L;wave_profile.depth],1.0);
                     depth9=interp1(result(1,:),result(2,:),wave_profile.L);  
                     depth9(1:T)=NaN;depth9(cliffjun_loc:end)=NaN;  
                     a=diff(depth9);b=diff(a);
                     [pks,locs] = findpeaks(b,0.01,1,0.2);
                     
                     if ~isempty(locs);
                         toe=locs(1)+1;
                         toe_onshore=toe;
                     
                         
                     end
                     
                     
                     end
                     
                  end

                    
                end
 
                end 
            end
        end
        
        
     end

     
%% Pick out maxi ensure it correctly corresponds to the maximum location 
 [pks,locs] = findpeaks(wave_profile.depth,[],1,0.01);
locs(locs<maxi-4)=[];
locs(locs>maxi+4)=[];

if ~isempty(locs);
    d=abs(locs-maxi);
    [~,d]=min(d);
    maxi=locs(d);
end

% if maxi onshore exists, its likely that the relevant mximum comes before
% a potential cliff
if exist('mod_toe','var') && exist('maxi_onshore','var');
    cliffjun_loc=maxi_onshore;
end


%Fix up where the cliffjun really isnt a slope break 

if strcmp(ESI,'2A') && wave_profile.depth(maxi)>10 &&  wave_profile.depth(cliffjun_loc)>9 && cliffjun_loc>30 && maxi-cliffjun_loc<8
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],1.0);
    depth11=interp1(result(1,:),result(2,:),wave_profile.L);  
    
    
    fixx=movingslope(depth11,5);
    fixx=movingslope(fixx,5);fixx=fixx*-1;
    [pks,locs] = findpeaks(fixx,[],1,0.01);
    if ~ismember(cliffjun_loc,locs)
        cliffjun_loc=maxi;
    end
        
end


if strcmp(ESI,'3A') && wave_profile.depth(maxi)>10 &&  wave_profile.depth(cliffjun_loc)>7.5 && cliffjun_loc>60
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.2);
    depth11=interp1(result(1,:),result(2,:),wave_profile.L);  
    
    fidx=atan(movingslope(depth11,5))*(180/pi);
    fidx=fidx(toe:maxi);
    
    A=find(fidx>12);
    
    if isempty(A);
    fixx=movingslope(depth11,5);
    fixx=movingslope(fixx,5);fixx=fixx*-1;
    [pks,locs] = findpeaks(fixx,[],1,0.01);
    locs(locs<toe)=[];
    if ~isempty(locs)
        cliffjun_loc=locs(1);
    else
        cliffjun_loc=maxi;
    end
    end
end




     
%% start prepping output
%maximum locaiton and elevation 
Out.maxi=maxi;
if ~isempty(Out.maxi)
Out.maxi_ele=wave_profile.depth(maxi);
else
Out.maxi_ele=[];
end

%If the calculated toe is lower than MHW or seaward of MHW location,
%replace with MHW parameters
if toe>0;
if wave_profile.depth(toe)<wave_profile.mhw && MHW_loc<toe;
    
    remp=wave_profile.depth;
    remp(1:toe)=NaN;
    temp=find(remp>wave_profile.mhw);
    if ~isempty(temp);
    temp=[temp(1)-1 temp(1)];
    if temp(1)<1
        I=1;
    else
    val=remp(temp);
    val=abs(val-wave_profile.mhw);
    [m I]=min(val);
    end
    
    MHW=wave_profile.mhw;
    MHW_loc=temp(I);
    else
    MHW=wave_profile.mhw;;
    MHW_loc=find(~isnan(remp));   
    MHW_loc=MHW_loc(1);
    end
    toe=MHW_loc;
    
end
else
    toe=MHW_loc;
end

if toe==cliffjun_loc && isempty(Rey);
    toe=MHW_loc;
end


% Secondary catch to fix the toe in rocky envrionments 
if ~strcmp(ESI,'3A') && ~strcmp(ESI,'4') && ~strcmp(ESI,'6B') && toe>5;
    if exist('inflect2','var')
        if inflect_dist<15
       if min(wave_profile.depth(inflect1:inflect2))<7 
        tmp=wave_profile.depth;
        Dvec = movingslope(tmp,3);
        rmp=diff(Dvec);
        [pks,locs] = findpeaks(rmp,[],[],0.25);
        locs(locs>inflect2+5)=[];
        locs(locs>maxi)=[];
        if~isempty(locs);
           toe=locs(end);
        end
            
            
       end
        end
    end

    %if the bench is shoreward of the toe, then the bench information is
    %unecessary in runup calculations    
    if strcmp(ESI,'2A')
        if bench_loc>toe
            bench_loc=[];
            bermwidth=[];
            bench_toe_loc=[];
        end
    end
end
    
%% Fix if the toe elevation is very large 4/14/2019
%first normal toe location
if wave_profile.depth(toe) > 7
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],1.0);
    depth8=interp1(result(1,:),result(2,:),wave_profile.L);
    tmp=depth8(1:toe);
    %find a new toe inflection from here
    Dvec = diff(tmp);
    rmp=diff(Dvec);
    [pks,locs] = findpeaks(rmp,[],1,0.25);
    if ~isempty(locs)
    locs=locs+1;
    locs(locs>toe)=[];
    tmp=wave_profile.depth(locs);
    tmp=find(tmp>7);
    locs(tmp)=[];
    
    if ~isempty(locs)
        toe=locs(end);
    else
        toe=1;
    end
    end
end

%second toe onshore location
if toe_onshore>0
if wave_profile.depth(toe_onshore) > 7
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],1.0);
    depth8=interp1(result(1,:),result(2,:),wave_profile.L);
    tmp=depth8(1:toe_onshore);
    %find a new toe inflection from here
    Dvec = diff(tmp);
    rmp=diff(Dvec);
    [pks,locs] = findpeaks(rmp,[],1,0.25);
    if ~isempty(locs)
    locs=locs+1;
    locs(locs>toe_onshore)=[];
    tmp=wave_profile.depth(locs);
    tmp=find(tmp>7);
    locs(tmp)=[];
    if ~isempty(locs)
        toe_onshore=locs(1);
    else
        toe_onshore=1;
    end
    end
end
else
    toe_onshore=1;
end
%% Select the appropaiate overtopping point 
%two options, cliff junc location and maxi location

%first calculate the slope between the toe and the junction
%elect the toe closer onshore
if toe>=toe_onshore
    qual=toe;
else
    qual=toe_onshore;
end

if qual>1;
    point1=qual;
else
    point1=1;
end

%now calc the slopes
%slope from toe to slope break
slopejun=(wave_profile.depth(cliffjun_loc)-wave_profile.depth(point1))/(cliffjun_loc-point1);
%slope form toe to maximum
slopemax=(wave_profile.depth(maxi)-wave_profile.depth(point1))/(maxi-point1);
%slope between slope break and maximum
slopeint=(wave_profile.depth(maxi)-wave_profile.depth(cliffjun_loc))/(maxi-cliffjun_loc);

slopejun=atan(slopejun)*(180/pi);
slopemax=atan(slopemax)*(180/pi);
slopeint=atan(slopeint)*(180/pi);

%now have the slopes in degrees for comparison

%% Determine Overtopping Point (thrresholds determined through testing and 
%  calibrated by region in a cascading series of if then statements
if slopemax>slopejun && maxi-cliffjun_loc<80 && toe<100 && ~strcmp(ESI,'6B') || slopemax>slopejun+10 && slopemax>20 && maxi-cliffjun_loc<50 && toe<120 && strcmp(ESI,'3A');
    overtop_point=maxi;
else
    %here compare the slopes
    slopediff=slopejun-slopemax;
    
    if slopediff < 10 && maxi-cliffjun_loc<10 && wave_profile.depth(maxi)-wave_profile.depth(cliffjun_loc)<4 || slopediff < 15 && maxi-cliffjun_loc<10 && wave_profile.depth(maxi)-wave_profile.depth(cliffjun_loc)<4 || slopediff < 10 && wave_profile.depth(maxi)<7 && maxi-cliffjun_loc<90; %added second half 4/16/2019
        overtop_point=maxi;
    elseif strcmp(ESI,'3A') && slopemax>50 && slopeint>45;
        overtop_point=maxi;
    elseif strcmp(ESI,'3A') && slopemax<10 && slopejun<15 && slopeint<5 && maxi-cliffjun_loc>45;
        overtop_point=maxi;
    elseif strcmp(ESI,'3A')
       if slopediff < 10 && maxi-cliffjun_loc<10 || slopediff < 15 && maxi-cliffjun_loc<20 && wave_profile.depth(maxi)-wave_profile.depth(cliffjun_loc)<4 || slopediff < 10 && wave_profile.depth(maxi)<7 && maxi-cliffjun_loc<90; %added second half 4/16/2019
          overtop_point=maxi;
       elseif slopeint<15 && maxi-cliffjun_loc<30 && wave_profile.depth(maxi)>10 && wave_profile.depth(cliffjun_loc)>7 && wave_profile.depth(maxi)-wave_profile.depth(cliffjun_loc)<6;
          overtop_point=maxi;  
       elseif slopeint<10 && maxi-cliffjun_loc<45 && wave_profile.depth(maxi)>10 && wave_profile.depth(cliffjun_loc)>6 && wave_profile.depth(maxi)-wave_profile.depth(cliffjun_loc)<6;
          overtop_point=maxi;  
       else
          overtop_point=cliffjun_loc;
       end
    elseif strcmp(ESI,'2A')  
        [pks,locs] = findpeaks(wave_profile.depth,[],10,3);
        locs(locs>10)=[];
        if ~isempty(locs)
            if ismember(cliffjun_loc,locs);
               overtop_point=maxi; 
            else
               overtop_point=cliffjun_loc; 
            end
        else
            overtop_point=cliffjun_loc; 
        end
    elseif strcmp(ESI,'1A') && slopemax>70;
        overtop_point=maxi;    
    elseif strcmp(ESI,'5') && slopemax>50 && slopeint>50;
        overtop_point=maxi;   
    elseif strcmp(ESI,'1B') && slopemax>50 && slopeint>50;
        overtop_point=maxi;  
    elseif strcmp(ESI,'6B') && abs(slopediff)<7 && slopemax>20 && maxi-cliffjun_loc<20 && slopeint>slopejun;
        overtop_point=maxi;    
        
    else
        overtop_point=cliffjun_loc;
    end
    
    
end

%% Include more catches for special scenarios regarding overtopping points
if wave_profile.depth(maxi) < 15 && maxi-cliffjun_loc < 15 && wave_profile.depth(maxi)>wave_profile.depth(cliffjun_loc); %Maybe 10m and not 15m???
    overtop_point=maxi;
else
    %Keep what was determined before
end

if wave_profile.depth(cliffjun_loc)<3 && wave_profile.depth(maxi)>5
    overtop_point=maxi;
    
elseif wave_profile.depth(cliffjun_loc)<5 && wave_profile.depth(maxi)>5 && strcmp(ESI,'2A')
    overtop_point=maxi;
end

if strcmp(ESI,'6D') && wave_profile.depth(cliffjun_loc)<5 %Just to deal with rocks in the way
    overtop_point=maxi;
end

if isempty(Rey) && cliffjun_loc-toe<5 && wave_profile.depth(cliffjun_loc)-wave_profile.depth(toe)<2 && overtop_point==cliffjun_loc && maxi-cliffjun_loc<50
    overtop_point=maxi;
end

if strcmp(ESI,'3A') && overtop_point==maxi && wave_profile.depth(cliffjun_loc)>6 && slopeint<15; %Just to deal with rocks in the way
    AAA=movingslope(wave_profile.depth,3);AAA=atan(AAA)*(180/pi); AAA(cliffjun_loc:end)=[];AAA(1:toe)=[];
    BBB=(wave_profile.depth(cliffjun_loc)-wave_profile.depth(toe))/(cliffjun_loc-toe);BBB=atan(BBB)*(180/pi);
    
    AAA=find(AAA>45);
    if ~isempty(AAA)
       overtop_point=cliffjun_loc; 
    end
    
end



%% ALast piece to deal with odd overtopping points

if strcmp(ESI,'4') && overtop_point-toe>100 && maxi-cliffjun_loc<10 && toe<25 && overtop_point<length(depth4)
    X=[1 overtop_point];
    V=[wave_profile.depth(1) wave_profile.depth(overtop_point)];
    Xq=[1:overtop_point];
    Vq = interp1(X,V,Xq);
    Values=depth4(1:overtop_point)-Vq;
    raft=mean(Values);
    if raft<-4
       overtop_point=maxi_onshore; 
    end
elseif strcmp(ESI,'4') && overtop_point-toe>100 && maxi-cliffjun_loc<10 && toe<25 && overtop_point>length(depth4)
    overtop_point=maxi_onshore; 
end

if exist('maxi_onshore','var')
if strcmp(ESI,'3A') && overtop_point==maxi && overtop_point>50 && maxi_onshore>overtop_point && wave_profile.depth(maxi_onshore)>wave_profile.depth(overtop_point)
        overtop_point=maxi_onshore; 
elseif strcmp(ESI,'3A') && overtop_point==maxi && overtop_point<25 && wave_profile.depth(maxi_onshore)>wave_profile.depth(overtop_point) && wave_profile.depth(overtop_point)<3.7
        overtop_point=maxi_onshore; 
end
end


%% Finalize Outputs
%MHW
Out.MHW=wave_profile.depth(MHW_loc);
Out.MHW_loc=MHW_loc;
Out.MSL=wave_profile.msl;

if cliffjun_loc>maxi
    cliffjun_loc=maxi;
end

% Clifftop and CliffJun are the same, was leftover from previous versions
% for compatibility
Out.CliffTop_Ele=wave_profile.depth(cliffjun_loc);
Out.CliffTop_Loc=cliffjun_loc;
Out.CliffJun_Loc=cliffjun_loc;
Out.CliffJun_Ele=wave_profile.depth(cliffjun_loc);
Out.toe_loc=toe; 
Out.toe_ele=wave_profile.depth(toe);

%if the toe location is first in the profile, likely an error in the
%original DEM extent, so convert to MHW for 3A and 5 beach profiles and -1
%for all others as this will be used to help determine plunging cliffs
if Out.toe_loc==1 && wave_profile.depth(1)>MSL;
    if strcmp(ESI,'3A') || strcmp(ESI,'5')
        Out.toe_ele=MHW;
        Out.toe_loc=MHW_loc;
    else
        Out.toe_ele=MSL;
        Out.toe_loc=-1;
    end
end

%to account for errors in the toe_onshore variable
if exist('toe_onshore','var');
Out.toe_onshore=toe_onshore; 
end
if exist('maxi_onshore','var');
Out.maxi_onshore=maxi_onshore; 
end

%if in a plunging cliff environment, make toe_onshore also equal to -1 
if exist('toe_onshore','var');
if Out.toe_loc==-1;
   Out.toe_onshore=-1;
else 
    Out.toe_onshore=Out.toe_onshore;
end
end 

Out.berm_loc=bench_loc;
Out.berm_width=bermwidth;
if ~isempty(bench_loc) && bench_toe_loc==0;
    bench_toe_loc=1;
    
    if Out.berm_loc==0;
        Out.berm_loc=bench_toe_loc+1;
    end
end
Out.berm_toe=bench_toe_loc;

if exist('toez','var')
if ~isempty(toez)
    toey=unique(toez);
    
    if length(toey)==1;
        if abs(toe-toey)<5
        toe=toey;
        Out.toe_loc=toe; 
        Out.toe_ele=wave_profile.depth(toe);
        end
    end
    
end
end

%Final toe check for 2A environments
if strcmp(ESI,'2A') && ~isempty(Out.berm_loc) && Out.toe_loc>3;
    
    clope=(wave_profile.depth(Out.toe_loc)-wave_profile.depth(Out.berm_loc))/(Out.toe_loc-Out.berm_loc);;
    if clope < -0.02
        Out.toe_loc=Out.berm_toe;
        Out.toe_ele=wave_profile.depth(Out.toe_loc);
    end
end



%remaining profile information
Out.ESI=ESI;% ESI
Out.extra_toe=extra_toe;%if a potential exta toe candidate was flagged (for use in review)
Out.Rey=Rey;% Flag for a seawall
Out.depth=wave_profile.depth;%elevation profile
Out.xpf=wave_profile.xpf; % X coordinates of all points along profile
Out.ypf=wave_profile.ypf; % Y coordinates of all points along profile
Out.L=wave_profile.L; % Distance along profile in meters 
Out.overtop_point=overtop_point; %location of the determined overtopping point along the elevaiton profile 

end
