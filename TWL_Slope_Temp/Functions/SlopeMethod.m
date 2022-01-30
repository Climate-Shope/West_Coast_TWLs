function [Out]=SlopeMethod(wave_profile,test_con,Wavecon,cut,ESIID,Rey,maxtide,fslope);
warning off

%Input Variables 

%Wave_Profile = The Wave profile parameters including elevation and MHW
%values


%test_con = Maximum wave conditions at the profile to determine maximum
%runup extent 

%Wavecon = Used just to grab the depth of the wave measurments

%cut = if need to edit the first points along a profile, cut from 1 to
%point specified. 


%Discliamer: this is a much simplified version from previous versions in an
%effort to expedite the run time. 

%% Variables to consider and modify
window=5; %For the smooth slope calculations
view_range=15+1;%meters to look for a better infleciton point around the coarse resolution

%% Prefilter the profile before calculations

if length(wave_profile.depth(1,:))==1;
    wave_profile.depth=wave_profile.depth';
end

if ~isempty(cut)
wave_profile.depth(1:cut)=[];
wave_profile.xpf(1:cut)=[];
wave_profile.ypf(1:cut)=[];
wave_profile.L=1:length(wave_profile.depth);
treb=1;
end
%Prefilter the profile to remove any inital decrease in slope
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

%% Determine MHW location and elevation
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


%% Curtail the length of the profile to a max length of 400 m as mentioned by K. Serafin in personal communication
if length(wave_profile.depth)>400 && wave_profile.depth(400)>2;
    wave_profile.depth(401:end)=[];
    wave_profile.xpf(401:end)=[];
    wave_profile.ypf(401:end)=[];
    wave_profile.L=1:length(wave_profile.depth);
end



%% Now calculate the morphological parameters
%%
%% START WITH THE CALCULTION OF THE MOST SHOREWARD MAXIMUM 

if isempty(Rey)

%first_simplify the profile
result = DouglasPeucker([wave_profile.L;wave_profile.depth],2.5);
simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);

%Find the most shoreward maximum of the simplified profile
[pks,locs] = findpeaks(simp_prof,[],25,MHW);
if isempty(locs) && max(wave_profile.depth)<10;
result = DouglasPeucker([wave_profile.L;wave_profile.depth],1.0);
simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);
[pks,locs] = findpeaks(simp_prof,[],25,MHW);
end

if ~isempty(locs);
    maxi=locs(1);
else
    maxi=[];
end

%Zzoom into the relevant area using the base profile
if ~isempty(maxi);
qwert_prof=wave_profile.depth;
qwert_prof([1:(locs(1)-view_range), (locs(1)+view_range):end])=NaN;
[qwert_pks,qwert_locs] = findpeaks(qwert_prof,[],5,MHW);

%This will consider the points that are also just shoreward from the points
rem=find(qwert_locs>(maxi+5));
qwert_locs(rem)=[];qwert_pks(rem)=[];

if isempty(qwert_locs)
    maxi=maxi;
else
    %Want the greater of the seaward points
   [~,temp]=max(qwert_pks); 
    qwert_locs=qwert_locs(temp);
    
    maxi=qwert_locs;
end
end

else

result = DouglasPeucker([wave_profile.L;wave_profile.depth],2.5);

simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);

%Find the most shoreward maximum of the simplified profile
[pks,locs] = findpeaks(simp_prof,[],25,MHW);
if isempty(locs) && max(wave_profile.depth)<10;
result = DouglasPeucker([wave_profile.L;wave_profile.depth],1.0);
simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);
[pks,locs] = findpeaks(simp_prof,[],25,MHW);
end

if ~isempty(locs);
    maxi=locs(1);
else
    maxi=[];
end

%Zoom into the relevant area using the base profile
if ~isempty(maxi);
qwert_prof=wave_profile.depth;
qwert_prof([1:(locs(1)-view_range), (locs(1)+view_range):end])=NaN;
[qwert_pks,qwert_locs] = findpeaks(qwert_prof,[],5,MHW);

%This will consider the points that are also just shoreward from the points
rem=find(qwert_locs>(maxi+5));
qwert_locs(rem)=[];qwert_pks(rem)=[];

if isempty(qwert_locs)
    maxi=maxi;
else
    %Want the greater of the seaward points
   [~,temp]=max(qwert_pks); 
    qwert_locs=qwert_locs(temp);
    
    maxi=qwert_locs;
end

    maxi=maxi(1);

   
    %if there is a deteached seawall that is present, then set it
    %to look back and forth by about 10m to find max
    Rey=Rey(1,:);
    xcomp=(Rey(:,2)-wave_profile.xpf).^2;
    ycomp=(Rey(:,3)-wave_profile.ypf).^2;
    distance=sqrt(xcomp+ycomp);
    [~,I]=min(abs(distance));
    
    if I<maxi
    
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.3);
    simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);
    temp=simp_prof;
    temp2=nan(1,length(wave_profile.depth));
    try
    temp2(I-10:I+10)=temp(I-10:I+10);
    catch 
        if I-5<=0;
        temp2(1:I+10)=temp(1:I+10); 
        else
        
        temp2(I-5:I+10)=temp(I-5:I+10); 
        end
    end
    %If there is a detached seawall present, then need to limit search area
    %to the region close to the interseciton point 
    try
    [qwert_pks,qwert_locs] = findpeaks(temp2,[],5,MHW);
    catch
        maxi=[];
    end
    
    if ~isempty(qwert_locs);
        [~,K]=max(qwert_pks);
        maxi=qwert_locs(K);
    else
        maxi=I;
    end
    else
    
        maxi=maxi;
    end
else
   Rey=Rey(1,:);
    xcomp=(Rey(:,2)-wave_profile.xpf).^2;
    ycomp=(Rey(:,3)-wave_profile.ypf).^2;
    distance=sqrt(xcomp+ycomp);
    [~,I]=min(abs(distance)); 
     result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.3);
    simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);
    
    if I>300
        Rey=[];
        [pks,locs] = findpeaks(simp_prof,[],5,MHW);
        maxi=locs(1);
    else
    
    
    temp=simp_prof;
    temp2=nan(1,length(wave_profile.depth));
    try
    temp2(I-10:I+10)=temp(I-10:I+10);
    catch 
        if I-5<=0;
        temp2(1:I+10)=temp(1:I+10); 
        else
        
        temp2(I-5:I+10)=temp(I-5:I+10); 
        end
    end

    try

    [qwert_pks,qwert_locs] = findpeaks(temp2,[],5,MHW);
    catch
        maxi=[];
    end
    
    if ~isempty(qwert_locs);
        [~,K]=max(qwert_pks);
        maxi=qwert_locs(K);
    else
        maxi=I;
    end
    end
end
end   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build in a limit that the maximum can be from the MHW location
%For now, limit it to 400m as discussed with K. Serafin

if ~isempty(maxi);
    temp=maxi-wave_profile.mhw_Lpos;
    if temp>400;
        maxi=[];
    end
end



%%
%% STEP 2 CALCULATION OF THE TOE OF THE STRUCTURE/DUNE/ETC

%Create a dummy profile so that very nearshore toes can be found
CoastX=wave_profile.xpf(1);
CoastY=wave_profile.ypf(1);
CoastZ=wave_profile.depth(1);

OffshoreX=Wavecon.x0;
OffshoreY=Wavecon.y0;
OffshoreZ=Wavecon.z0-wave_profile.msl; %Referenced to MSL so subtract MSL to get NAVD88

%Now calculate the distance between the points
xcomp=(CoastX-OffshoreX)^2;
ycomp=(CoastY-OffshoreY)^2;
distance=round(sqrt(xcomp+ycomp));

%Now to interpolate between the two 
x=[1 distance];
y=[-1*OffshoreZ CoastZ];    
newx=1:distance;
int_bath=interp1(x,y,newx);
int_bath(end)=[];

distance=distance-1;

tmp=fliplr(newx(1:end-1));
tmp=wave_profile.L(1)-tmp;
tmp=[tmp wave_profile.L];
int_bath=[int_bath wave_profile.depth];

result = DouglasPeucker([tmp;int_bath],1);
simp_prof=interp1(result(1,:),result(2,:),tmp);

%now if maxi exists, the toe should be seaward by definition, so turn
%everything greater than maxi into NaN values
if ~isempty(maxi)
    simp_prof(maxi+distance+1:end)=NaN;
end

%do not want to cut to MHW in case the toe is less than MHW (which it
%shouldn't be for sandy coasts, but just in case for the rocky coastlines)

Dvec_simp = movingslope(simp_prof,window); 
Dvec2_simp = movingslope(Dvec_simp,window);

try
[pks_d2,locs_d2] = findpeaks(Dvec2_simp,0.0001,3,0.01);
catch
    locs_d2=[];
end

if ~isempty(locs_d2);
    locs_d2(locs_d2<(MHW_loc+distance))=[]; %At this point, can say that if it is less than MHW, dont want it
end

Dvec_temp = movingslope(simp_prof,2); %To make sure the slope measurement is highly localized
 %This looks at the local area and determines the slope of the potential
 %toe locations. If greater than 45 degrees, it is likely just a slope
 %change along a cliff or armoring structure. 


tmp_slope=Dvec_temp(locs_d2);

locs_d2(tmp_slope>0.8391)=[]; 
locs_d2(locs_d2>200+distance)=[];

uyp=locs_d2;

locs_d2(locs_d2<MHW_loc+distance)=[];

if isempty(locs_d2) && ~isempty(uyp);
    locs_d2=uyp(1);
end

     if ~isempty(locs_d2)
    temp=Dvec_temp(locs_d2);
    temp=atan(temp).*(180/pi);
    temp=find(temp>45);
    locs_d2(temp)=[];
    end

%This results in the coarser simplification, so iterate Stockdon if
%there are multiple toe locations. Think of this as the first pass,
%however, if there is only one toe, then there is only one toe.

%%%%%%%%%%%%%%%%%%%%%% ABuid in Another Fix %%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ESIID,'3A') && numel(locs_d2)==1 && locs_d2-distance <=1;
    locs_d2=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(locs_d2);
    
elseif numel(locs_d2)==1;

else
    
%If there are multiple potential locations, use stockdon method as a first pass
%First load in the correct depth of the wave measurments
    
%%

%start by calculating the slopes

locs_d2=locs_d2-distance;

 comp1=(wave_profile.depth(locs_d2)-wave_profile.depth(MHW_loc));
    comp2=locs_d2-MHW_loc;
    slopes=comp1./comp2;
    slopes=atan(slopes);
slopes((locs_d2-distance)<0)=atan(fslope(:,3));
%remove the slopes that are greater than 34 degrees (where stockdon et al, 2006 runup method would
%not work, this indicating that the TAW approach must be used).
%Any greater than 34 degrees is likely to be an erronous point in any case
    temp=find((slopes.*(180/pi))>34);
    slopes(temp)=[];
    locs_d2(temp)=[];

[R,out1,out2] = Runup(test_con.Hs,test_con.T,slopes,13); R= R+wave_profile.msl+maxtide; % Add runup to the still water level
    ver=[];
 if ~isempty(R);
    for mm=1:length(R);
    temp=find(R(mm)>wave_profile.depth(locs_d2(mm)));
    if ~isempty(temp)
        temp=mm;
    else
        temp=[];
    end
    ver=[ver,temp];
    end
    
    if  isempty(ver);
        locs_d2=locs_d2(1); %If it does not reach the elevation fo the toe, then just assume that the toe is the first one in series 
    else
        
        if numel(ver)==numel(R);
            [~,temp]=max(R); %If multiple, want the one that has the greatest slope (i.e. runup)
            locs_d2=locs_d2(temp);
        else
            locs_d2=locs_d2(ver);
        end
        

    end
 end
end
%At this point, there should be a limited (if not empty) number of
%locations where the infleciton points from the profile are located 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if it is empty, then need to double check using a less simplified version
%of the profile. 

if isempty(locs_d2);
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.2);
    simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);

    %now if maxi exists, the toe should be seaward by definition, so turn
    %everything greater than maxi into NaN values
    if ~isempty(maxi)
        simp_prof(maxi+1:end)=NaN;
    end

    Dvec_simp = movingslope(simp_prof,window); 
    Dvec2_simp = movingslope(Dvec_simp,window); 
    [pks_d2,locs_d2] = findpeaks(Dvec2_simp,0.001,5,[]);
    %ideally, the toe should not be less than the MHW mark
    pks_d2(locs_d2<MHW_loc)=[];
    locs_d2(locs_d2<MHW_loc)=[];

    Dvec_temp = movingslope(simp_prof,2); %To make sure the slope measurement is highly localized
    if ~isempty(locs_d2)
    temp=Dvec_temp(locs_d2);
    temp=atan(temp).*(180/pi);
    %if the slope is too great at the potential toe, remove the location
    temp=find(temp>45);
    locs_d2(temp)=[];
    pks_d2(temp)=[];
    end
    
    
    % remove location where the pks_d2 is negative in this instance
    pks_d2(pks_d2>0)=[];
    locs_d2(pks_d2>0)=[];
    
    %now remove the locaitons where the toe is greater than say 6 m higher
    %than the MHW (determined by regional testing)
    if ~isempty(locs_d2)
    tmp_el=wave_profile.depth(locs_d2)-MHW;
    locs_d2(tmp_el>6)=[];
    end
       
    %Still no Toe?
    
    if isempty(locs_d2);
       locs_d2=MHW_loc; %If it is still empty, assume that there isnt a toe, so call it MHW
    else %check the locations versus calculated runup 

%If there are multiple potential locations, use stockdon as a first pass
    
%%
%start by calculating the slopes
    comp1=(wave_profile.depth(locs_d2)-wave_profile.depth(MHW_loc));
    comp2=locs_d2-MHW_loc;
    slopes=comp1./comp2;
    slopes=atan(slopes);
    slopes((locs_d2-distance)<0)=atan(fslope(:,3));
%remove the slopes that are greater than 34 degrees (where stockdon et al (2006) runup method would
%not work, this indicating that the TAW approach must be used).
%Any greater than 34 degrees is likely to be an erronous point in any case
    temp=find((slopes.*(180/pi))>34);
    slopes(temp)=[];
    locs_d2(temp)=[];

    %once again do the runup test
    
 [R,out1,out2] = Runup(test_con.Hs,test_con.T,slopes,13); R= R+wave_profile.msl+maxtide;
    ver=[];
 if ~isempty(R);
    for mm=1:length(R);
    temp=find(R(mm)>wave_profile.depth(locs_d2(mm)));
    if ~isempty(temp)
        temp=mm;
    else
        temp=[];
    end
    ver=[ver,temp];
    end
    
    if  isempty(ver);
        locs_d2=locs_d2(1); %If it does not reach the elevation fo the toe, then just assume that the toe is the first one in series 
    else
        locs_d2=locs_d2(ver);
    end
 end
 
if isempty(ver)&& ~isempty(locs_d2);
        locs_d2=locs_d2(1); %If it does not reach the elevation fo the toe, then just assume that the toe is the first one in series 
    elseif isempty(ver)&& isempty(locs_d2);
        locs_d2=MHW_loc;

    else
        locs_d2=locs_d2(end); %If multiple points are inundated, want to grab the last one in series
    end 
    
    end
    view_range=20+1;%meters
    
else
    view_range=15+1;%meters
    
    
%Create a dummy profile so that very nearshore toes can be found
CoastX=wave_profile.xpf(1);
CoastY=wave_profile.ypf(1);
CoastZ=wave_profile.depth(1);

OffshoreX=Wavecon.x0;
OffshoreY=Wavecon.y0;
OffshoreZ=Wavecon.z0-wave_profile.msl; %Referenced to MSL so subtract MSL to get NAVD88

%Now calculate the distance between the points
xcomp=(CoastX-OffshoreX)^2;
ycomp=(CoastY-OffshoreY)^2;
distance=round(sqrt(xcomp+ycomp));

%Now to interpolate between the two 
x=[1 distance];
y=[-1*OffshoreZ CoastZ];    
newx=1:distance;
int_bath=interp1(x,y,newx);
int_bath(end)=[];

tmp=fliplr(newx(1:end-1));
tmp=wave_profile.L(1)-tmp;
tmp=[tmp wave_profile.L];
int_bath=[int_bath wave_profile.depth];
      
   distance=distance-1; 

   if  locs_d2(1)<distance;
       locs_d2=locs_d2+distance;
   end
    jj=[0.1];window=3; %Changes for the smaller jj constraint 
    qwert = DouglasPeucker([tmp;int_bath],jj);
    qwert_prof=interp1(qwert(1,:),qwert(2,:),tmp);
    temp=NaN(1,length(qwert_prof));
    range=[];
    for ii=1:length(locs_d2);
        if ii==length(locs_d2);
        r=[(locs_d2(ii)-view_range):(locs_d2(ii)+view_range)];
        else
        r=[(locs_d2(ii)-view_range):(locs_d2(ii)+view_range)];
        end
        range=[range,r];
    end
    range=unique(range);
    range(range<1)=[];
    if ~isempty(maxi+distance);
    range(range>(maxi+distance))=[];
    end
    range(range>length(qwert_prof))=[];
    temp(range)=qwert_prof(range);
    qwert_prof=temp;
    Dvec_qwert = movingslope(qwert_prof,window); 
    Dvec2_qwert = movingslope(Dvec_qwert,window); 
    

    [qwert_pks_d2,qwert_locs_d2] = findpeaks(Dvec2_qwert,0.001,3,0.001);
    if isempty(qwert_locs_d2);
        [qwert_pks_d2,qwert_locs_d2]=max(Dvec2_qwert);
    end
    
    qwert_locs_d2=[qwert_locs_d2,locs_d2];
    qwert_locs_d2=sort(qwert_locs_d2); %Puts it in the right order
    
%%%%%%%%%%%%%%%%%%%%%% Attempting to Rectify %%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(ESIID,'3A');
    qwert_locs_d2=qwert_locs_d2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Dvec_temp = movingslope(qwert_prof,2);
if ~isempty(qwert_locs_d2)
    temp=Dvec_temp(qwert_locs_d2);
    temp=atan(temp).*(180/pi);
    temp=find(temp>45);
    
    
    if ~isempty(temp) %Look at the preceding point to see if it better fits the bill
    qwert_locs_d2(temp)=qwert_locs_d2(temp)-1;
    temp=Dvec_temp(qwert_locs_d2);
    temp=atan(temp).*(180/pi);
    temp=find(temp>45);
    end
    
    
    qwert_locs_d2(temp)=[];
end
    
 if ~isempty(qwert_locs_d2)
   
    quan=find(int_bath(qwert_locs_d2)>(10*MHW));% arbitrary cutoff, but keeps the toe from being unrealistically high
        qwert_locs_d2(quan)=[];   
    
     if isempty(qwert_locs_d2);
            qwert_locs_d2=MHW_loc+distance;
     end
     
     if numel(qwert_locs_d2)==1 && qwert_locs_d2==MHW_loc+distance;
        temp_loc=qwert_locs_d2;
    else     %calculate runup for these locations
  
    % Calculate slopes for stockdon runup
    comp1=(int_bath(qwert_locs_d2)-int_bath(MHW_loc+distance));
    comp2=qwert_locs_d2-MHW_loc+distance;
    slopes=comp1./comp2;
    slopes=atan(slopes);
     
    slopes((qwert_locs_d2-distance)<0)=atan(fslope(:,3));
        
    
    
    %remove the slopes that are greater than 34 degrees (where stockdon et al. (2006) runup method would
    %not work, this indicating that the TAW approach must be used).
    %Any greater than 34 degrees is liekly to be an erronous point in any case
    temp=find((slopes.*(180/pi))>34);
    slopes(temp)=[];  
    
    %once again do the runup test
    [R,out1,out2] = Runup(test_con.Hs,test_con.T,slopes,13); R= R+wave_profile.msl+maxtide;
    ver=[]; 
     
    for mm=1:length(R);
    temp=find(R(mm)>int_bath(qwert_locs_d2(mm)));
     if ~isempty(temp)
        temp=mm;
     else
        temp=[];
     end
     ver=[ver,temp];
    end
    
    if isempty(ver)&& ~isempty(qwert_locs_d2);
        temp_loc=qwert_locs_d2(1); %If it does not reach the elevation fo the toe, then just assume that the toe is the first one in series 
    elseif isempty(ver)&& isempty(qwert_locs_d2);
        temp_loc=MHW_loc+distance;
    
    else
        qwert_locs_d2=qwert_locs_d2(ver);
        temp_loc=qwert_locs_d2(end); %If multiple points are inundated, want to grab the last one in series
    end 
       
    end    
     
else
    temp_loc=MHW_loc+distance;
 end     

 locs_d2=temp_loc;
 locs_d2=locs_d2-distance;
 if locs_d2==0
     locs_d2=-1;
 end
end 
 
if isempty(locs_d2); %here would be an almost vertical slope ot of the water to a maximum
       toe_loc=MHW_loc;
elseif numel(locs_d2)==1;
    %at this point, there should only be one toe location identified
       toe_loc=locs_d2;
else
    %if for any reason, more than one has been returned just grab the last
    %one as the runup level should still be higher
       toe_loc=locs_d2(end);
end  
 

%% If there is a detached seawall then do the same thing, look around the seawall location, and find the shoeward toe
    if ~isempty(Rey);
        if maxi>toe_loc+100;
            toe_loc=toe_loc;
        else
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.3);
    temp=interp1(result(1,:),result(2,:),wave_profile.L);

    temp2=nan(1,length(wave_profile.depth));
    try
    temp2(maxi-10:maxi+10)=temp(maxi-10:maxi+10);
    catch
        if maxi-5<=0;
            temp2(1:maxi+10)=temp(1:maxi+10);  
        else
            temp2(maxi-5:maxi+10)=temp(maxi-5:maxi+10);  
        end
    end
    
    Dvec_simp = movingslope(temp2,3); 
    Dvec2_simp = movingslope(Dvec_simp,3); 
    [pks_d2,locs_d2] = findpeaks(Dvec2_simp,0.001,5,[]);
    try
        locs_d2=locs_d2(1);
        toe_loc=locs_d2;
    catch
        toe_loc=maxi-1;
    end
        end
    end

%% Now Add in the same common sense parameter where the toe should not be above 6m
try 
    if (wave_profile.depth(toe_loc)-MHW)>6;
    toe_loc=MHW_loc;
    end
catch

end


%% Step 3: Now determine if there is a berm seaward of the toe location 
%the next step will be to dermine if need to find inflection point 
%landward of the toe location 

%Finding the berm in this case should be easy if MHW and the TOE exist, if
%not, things get a little more complicated 

%first edit the profile in a simple way
%Start by simplifying the profile again, but this time not as coarse
if toe_loc~=MHW_loc;
result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.08);
simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);

%now if maxi exists, the toe should be seaward by definition, so turn
%everything greater than maxi into NaN values
if ~isempty(maxi)
    simp_prof(maxi+1:end)=NaN;
end

if wave_profile.mhw_Lpos < toe_loc;

    if wave_profile.mhw_Lpos > 0;
    simp_prof(1:wave_profile.mhw_Lpos)=NaN;
    end

    simp_prof(toe_loc+1:end)=NaN;
end    

if wave_profile.mhw_Lpos > toe_loc;

    if wave_profile.mhw_Lpos > 0;
    simp_prof(wave_profile.mhw_Lpos+1:end)=NaN;
    end

    simp_prof(1:toe_loc)=NaN;
end   




%do not want to cut to MHW in case the toe is less than MHW (which it
%shouldn't be for sandy coasts, but just in case for the rocky coastlines)
if numel(find(~isnan(simp_prof)))<15;
    window=3;
end
Dvec_simp = movingslope(simp_prof,window); 
Dvec2_simp = movingslope(Dvec_simp,window);

if numel(find(~isnan(simp_prof)))<8;
Dvec_simp = diff(simp_prof); 
Dvec2_simp = diff(Dvec_simp);
temp=Dvec2_simp*-1;
try
[pks_d2,locs_d2] = findpeaks(temp,0.0001,3,0.001);
locs_d2=locs_d2+1;
catch
    disp('no berm');
 locs_d2=[];
end

else
%okay, the simplest way after cutting between the MHW and toe location,
%should be the greatest peak there. 
temp=Dvec2_simp*-1;
try
[pks_d2,locs_d2] = findpeaks(temp,0.0001,3,0.001);
catch
    disp('no berm');
 locs_d2=[];
end
end


%add in a clause to look for a slope greater than 32 degrees 
if ~isempty(locs_d2);
    if numel(locs_d2)==1;
       berm_loc= locs_d2;
    else
       [~,temp]=max(pks_d2);
       berm_loc= locs_d2(temp);
    end
    
    
    compx=berm_loc-wave_profile.mhw_Lpos;
    compy=wave_profile.depth(berm_loc)-wave_profile.depth(MHW_loc);
    slope=compy/compx;
    slope=atan(slope);
    if (slope*(180/pi))<32 && ~strcmp(ESIID,'2B');
        berm_loc=0;
    end
        
    
    
    
    
else
    berm_loc=0;
end

else
    berm_loc=0;
end


if berm_loc~=0;
    
    if berm_loc < toe_loc
        berm_width=toe_loc-berm_loc;
    else
        Dvec = movingslope(wave_profile.depth,window);
        Dvec2 = movingslope(Dvec,window);
        Dvec2(1:berm_loc)=NaN; %Only concerned with whats coming next in order
        [pks_d2,locs_d2] = findpeaks(Dvec2,0.001,3,0.0001);
        locs_d2=locs_d2(1); %it should be the next one in order
        berm_width=locs_d2-berm_loc;
    end
else
    berm_width=0;
end


%Plug in a common sense value value to determine the berm, esentially, it
%should be floating around the toe, not much greater. 

if berm_loc ~= 0; 
    if toe_loc<1;
       toe_ele=MSL;
    else
       toe_ele= wave_profile.depth(toe_loc);
    end
    if wave_profile.depth(berm_loc) > toe_ele+6;
        berm_loc = 0;
    end
end

if toe_loc<1;
    berm_loc = 0;
end

%% Now determine if there is a cliff edge for overtopping
% determing a cliff edge should be as simple as determining the presence of
% a berm, just at a larger scale a
% Must occur after the TOE, if present.

result = DouglasPeucker([wave_profile.L;wave_profile.depth],4);
simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);

%Remove everything before the toe
simp_prof(1:toe_loc)=NaN;

Dvec_simp = movingslope(simp_prof,window); 
Dvec2_simp = movingslope(Dvec_simp,window);

temp=Dvec2_simp*-1;
try
    [pks_d2,locs_d2] = findpeaks(temp,0.0001,3,0.01);
catch
    locs_d2=[];
end

%for now, assume that if there are multiple points, want the one closer
%to shore. 

if ~isempty(locs_d2)
    locs_d2=locs_d2(1);
    
    CliffTop_Loc=locs_d2;
    CliffTop_Ele=wave_profile.depth(locs_d2);
    
else
    CliffTop_Loc=0;
    CliffTop_Ele=0;
end

%Actually grab a point a bit more seaward form the true top of the
%cliff to define TAW, will call this the Cliff Junction 

%% Cliff Junction calculation 
if CliffTop_Loc ~= 0;
result = DouglasPeucker([wave_profile.L;wave_profile.depth],2);
simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);

%Remove everything before the toe

dummy=repmat(50,1,10);
simp_prof=[dummy, simp_prof];
Dvec_simp = movingslope(simp_prof,3); 
Dvec2_simp = movingslope(Dvec_simp,3);

temp=Dvec2_simp*-1;
try
    [pks_d2,locs_d2] = findpeaks(temp,0.0001,3,0.01);
catch
    locs_d2=[];
end
locs_d2=locs_d2-length(dummy);

pks_d2(locs_d2<1)=[];
locs_d2(locs_d2<1)=[];

pks_d2(locs_d2<toe_loc)=[];
locs_d2(locs_d2<toe_loc)=[];

pks_d2(locs_d2>=CliffTop_Loc)=[];
locs_d2(locs_d2>=CliffTop_Loc)=[];

%for now, assume that if there are multiple points, we want the one closer
%to shore. 

if ~isempty(locs_d2)
    [~,I]=max(pks_d2);
    locs_d2=locs_d2(I);
    
    CliffJun_Loc=locs_d2;
    CliffJun_Ele=wave_profile.depth(locs_d2);
    
else
    CliffJun_Loc=CliffTop_Loc;
    CliffJun_Ele=CliffTop_Ele;
end
else
    CliffJun_Loc=CliffTop_Loc;
    CliffJun_Ele=CliffTop_Ele;  
end
%% if this is a detached seawall location, want to repllace the clifftop location with the crest of the seawall 
if ~isempty(Rey);
result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.3);
    simp_prof=interp1(result(1,:),result(2,:),wave_profile.L);
    temp=simp_prof;
    temp2=nan(1,length(wave_profile.depth));
    try
    temp2(maxi-10:maxi+10)=temp(maxi-10:maxi+10);
    catch 
        if maxi-5<=0;
            temp2(1:maxi+10)=temp(1:maxi+10);  
        else
            temp2(maxi-5:maxi+10)=temp(maxi-5:maxi+10); 
        end
    end
    
    Dvec_simp = movingslope(temp2,window); 
    Dvec2_simp = movingslope(Dvec_simp,window);
    temp=Dvec2_simp*-1;
    try
    [pks_d2,locs_d2] = findpeaks(temp,0.0001,3,0.01);
    catch
    locs_d2=[];
    end
    
    if ~isempty(locs_d2)
    locs_d2=locs_d2(1);
    
    CliffTop_Loc=locs_d2;
    CliffTop_Ele=wave_profile.depth(locs_d2);
    
    else %if it cannot be found, then just try to grab where the intersection point would be
    CliffTop_Loc=maxi;
    CliffTop_Ele=wave_profile.depth(maxi);
    end
    
    if CliffTop_Loc<toe_loc;
        %if the crest of the seawall is seaward of the toe, that means that
        %the toe needs to be shifted. 
        
    result = DouglasPeucker([wave_profile.L;wave_profile.depth],0.3);
    temp=interp1(result(1,:),result(2,:),wave_profile.L);

    temp2=nan(1,length(wave_profile.depth));
    try
    temp2(maxi-10:maxi+10)=temp(maxi-10:maxi+10);
    catch
        if maxi-5<=0;
           temp2(1:maxi+10)=temp(1:maxi+10); 
        else
           temp2(maxi-5:maxi+10)=temp(maxi-5:maxi+10);
        end
    end
    
    Dvec_simp = movingslope(temp2,3); 
    Dvec2_simp = movingslope(Dvec_simp,3); 
    [pks_d2,locs_d2] = findpeaks(Dvec2_simp,0.001,5,[]);
    zemp=find(locs_d2>CliffTop_Loc);
    locs_d2(zemp)=[];
    
     try
        locs_d2=locs_d2(1);
        toe_loc=locs_d2;
    catch
        toe_loc=maxi-1;
    end
    

    end
    
end
    

%% Rectify the location of the MHW mark in relation to the toe. 

if toe_loc<1;
    toe_loc2=toe_loc;
    toe_loc=1;
end

if wave_profile.depth(toe_loc)<wave_profile.mhw && MHW_loc<toe_loc;
    
    remp=wave_profile.depth;
    remp(1:toe_loc)=NaN;
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
    toe_loc=MHW_loc;
    
end


Out.MHW=wave_profile.depth(MHW_loc);
Out.MHW_loc=MHW_loc+addval;
Out.MSL=wave_profile.msl;

Out.CliffTop_Ele=CliffTop_Ele;
if CliffTop_Loc>0;
    Out.CliffTop_Loc=CliffTop_Loc+addval; 
else
    Out.CliffTop_Loc=CliffTop_Loc;
end

if CliffJun_Loc>0;
    Out.CliffJun_Loc=CliffJun_Loc+addval; 
else
    Out.CliffJun_Loc=CliffJun_Loc;
end


if CliffTop_Loc<CliffJun_Loc; %Basically if there is a seawall detected that resets the clifftop data

Out.CliffJun_Ele=CliffTop_Ele;  
Out.CliffJun_Loc=CliffTop_Loc;  

else
Out.CliffJun_Ele=CliffJun_Ele;  
Out.CliffJun_Loc=CliffJun_Loc;   
end

%Out.toe_loc=toe_loc;
if exist('toe_loc2','var');
if toe_loc2<1;
    Out.toe_ele=MSL;
else
    
    Out.toe_ele=wave_profile.depth(toe_loc2);
    
end

if toe_loc2>0;
    toe_loc2=toe_loc2+addval;
end

Out.toe_loc=toe_loc2;
else
    if toe_loc<1;
    Out.toe_ele=MSL;
    else
    Out.toe_ele=wave_profile.depth(toe_loc);
    end
    
    if toe_loc>0;
        toe_loc=toe_loc+addval;
    end
    
    Out.toe_loc=toe_loc;
end



if Out.toe_loc==1 && wave_profile.depth(1)>MSL;
    Out.toe_ele=MSL;
    Out.toe_loc=-1;
end




Out.maxi=maxi;
if ~isempty(Out.maxi)
Out.maxi_ele=wave_profile.depth(maxi);
else
Out.maxi_ele=[];
end

if ~isempty(maxi);
   maxi=maxi+addval;
end
Out.maxi=maxi;

if berm_loc>0;
   berm_loc=berm_loc+addval;
end
    
Out.berm_loc=berm_loc;
Out.berm_width=berm_width;
Out.ESI=ESIID;
    
