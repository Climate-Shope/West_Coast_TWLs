%Code to remove islands and rocks from the start of the Major 100-m spaced
%profiles

warning off all 
%% CONFIG 
NAME_GRD_transects  = 'Douglas'; 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');
%% 
VISIBILITY = 'on'

Above_Limit=30; %Horizontal distance above msl that must be exceeded for the profile section to be considered a non island or rock

%% 

dirin       = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_profiles_wave_model']; 
dirout      = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\interp_profiles_wave_model_edited']; if ~exist([dirout],'dir'), mkdir(dirout); end


load(['F:\West_Coast_TWL_Hazards\01_Data\Transects',filesep, NAME_GRD_transects,'_Transects_v3.mat']);
eval(['shoreline=',NAME_GRD_transects,'_Transects_v3;']);
X0=shoreline.X(:,2);
Y0=shoreline.Y(:,2);


load('Transects2.mat');
Ntr=1:length(Transects2.xpf(:,1));

%%

for ii=Ntr;
    disp([num2str(ii) '/' num2str(length(Ntr))])
    OBJECTID=ii;
    name = [dirin filesep 'profile_', dec2base(OBJECTID,10,4) ,'.mat'];
    if ~exist(name,'file');
        continue
    end
    
    load(name);

    xx=X0(OBJECTID);
    yy=Y0(OBJECTID);
    
    %for the profile, identify where the profiles are greater
    %than MSL for a short duration 
    
    mhw=wave_profile.mhw;
    msl=wave_profile.msl;
    fmw=msl;
    depth=wave_profile.depth;
    x=1:length(depth);
    
    L1=[x;depth];
    L2=[[1,length(depth)];[fmw,fmw]];
    P = InterX(L1,L2);
    P=round(P(1,:));
    
    if depth(1)>fmw
       P=[1,P];
    end
    
    temp=diff(P);
    if isempty(temp); %This means there is only one or no intersection so can skip
         depth2=depth;
    else
    
    if mod(length(P),2)==1; %add to the end if there isnt another crossing with the line further inland
        P=[P,length(P)+75]; %add a lot more to ensure that it does meet the correct conditions
    end
    
    
    %now loop through P and cut out the depth less than 75 until
    %reaching the limit
    depth2=depth;
    for jj=1:2:length(temp);
        if temp(jj)>Above_Limit;
                idx=P(jj)-1;
                if idx>10
                depth2(1:idx)=NaN;
                end
                break
            
        else 
                idx=P(jj+1);
                depth2(1:idx)=NaN;
         
        end
    end
    end 
%now find the nans and remove them from all measurments

temp=isnan(depth2);
temp=find(temp);

wave_profile.depth(temp)=[];
wave_profile.xpf(temp)=[];
wave_profile.ypf(temp)=[];
wave_profile.L(temp)=[];


%Okay, want to rerun the intersection code, and then grab the point closest
%to the shoreline

    depth=wave_profile.depth;
    x=1:length(depth);
    
    L1=[x;depth];
    L2=[[1,length(depth)];[fmw,fmw]];
    P = InterX(L1,L2);
    P=round(P(1,:));
    
    if depth(1)>fmw
       P=[1,P];
    end
    
    if ~isempty(P) && length(P)>1;
    xloc=wave_profile.xpf(P);
    yloc=wave_profile.ypf(P);
    
    xcomp=(xloc-xx).^2;
    ycomp=(yloc-yy).^2;
    distance=sqrt(xcomp+ycomp);
    [~,I]=min(distance);
    
    wave_profile.depth(1:P(I))=[];
    wave_profile.xpf(1:P(I))=[];
    wave_profile.ypf(1:P(I))=[];
    wave_profile.L(1:P(I))=[];
    elseif length(P)==1;
        wave_profile.depth(1:P)=[];
        wave_profile.xpf(1:P)=[];
        wave_profile.ypf(1:P)=[];
        wave_profile.L(1:P)=[];
    end





mhw=wave_profile.mhw;

 temp=find(wave_profile.depth>mhw);
    if ~isempty(temp);
    temp=[temp(1)-1 temp(1)];
    if temp(1)<1
        I=1;
    else
    val=wave_profile.depth(temp);
    val=abs(val-mhw);
    [m I]=min(val);
    end
    
    wave_profile.mhw=mhw;
    wave_profile.mhw_Lpos=temp(I);
    else
     wave_profile.mhw=mhw;
    wave_profile.mhw_Lpos=wave_profile.L(end);   
    end
    
 msl=wave_profile.msl;    
    temp=find(wave_profile.depth>msl);
    if ~isempty(temp);
    temp=[temp(1)-1 temp(1)];
    if temp(1)<1
        I=1;
    else
    val=wave_profile.depth(temp);
    val=abs(val-msl);
    [m I]=min(val);   
    end
    wave_profile.msl=msl;
    wave_profile.msl_Lpos=temp(I);
    else
        wave_profile.msl=msl;
        wave_profile.msl_Lpos=wave_profile.L(end);
    end


    fout = [dirout, filesep, 'profile_', dec2base(OBJECTID,10,4) ,'.mat']; 
    m=matfile(sprintf('%s', fout),'writable',true); 
    m.wave_profile = wave_profile; 
end
    
 winopen(dirout)
