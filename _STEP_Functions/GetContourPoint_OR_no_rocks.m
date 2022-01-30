function I = GetContourPoint_OR_no_rocks(profile_bathy,contour_depth,gap_tolerance,prof_res);


%GetContourPoint tracks the points along a bathy profile and pulls the index
%of the poiont that most closesly matches a desirect contour depth and if
%there are multiple instances, the closest point to the depth offshore. 
%The function is also designed to handle rocks and small islands, instead
%grabbing the index of the deepest point between the rock/island and the
%mainland. This point will be used to back calculate 15 m water depth
%conditions via a conservation of energy relationship. 

%profile_bathy = 1-dimensional input of all of the depths along the profile
%depths should be positive

%contour_depth = contour depth of interest for exaple 15 m. Value should be
%positive

%gap_tolerance = the maximum allowable distance between points above sea
%level for these points to be considered part of the same land block. The
%goal of the tolerance is to filter any small spits or creeks that would
%result in erroneous values. For example, a value of 10 means that any above
%water points separated by 10 meters or less will be considered a
%continuous block. 

%prof_res = resolution of the bathy profile (1m, 2m, etc.) 


dep_prof=profile_bathy;
h1=contour_depth;
dif_tol=gap_tolerance;
dif_tol=dif_tol/prof_res;


%function to identify the 15m water depth point 

%tolerance parameters;

%IF THE WHOLE PROFILE IS UNDERWATER THEN CAN JUST TAKE THE MOST OFFSHORE
%15M CONTOUR

%IF THERE IS ANY LAND INVOLVED, MARCH FROM SHORE TO THE DEEPEST
%PART OF THE PROFILE BEFORE HITTING LAND. TAKES CARE OF ROCK AND BAYS. 
%THEN USE BACKCALC TO GET THE 15 M WAVE CONDITIONS 


%IF THE ENDPOINT OF the PROFILE IS LESS THAN 15 M, IT WILL
%GENERALLY MEAN IN A BAY OR SHALLOW WATER, so just grab the lowest point

%STEP 1 is the endpoint above water? 

%simplify, convert every point in depth profile that is less than 1 to -1

EP=(dep_prof(end)<0);

if EP==0 %if the end is not above sea level
    Above=find(dep_prof<0);
    
    if isempty(Above)
        if dep_prof(1)>h1 %if its deeper, likely mean on open coast
        sub_prof=dep_prof-h1;
        [Nah ind]=min(abs(sub_prof));
        cut=dep_prof(ind);    
        
        
        dep=find(dep_prof<=cut);
        I=dep(1);
        else % if not, likely in a bay and want the closest to 15m
            %first see if ther eis an intersection between the two 
            t=1:length(dep_prof);
            intx=repmat(h1,1,length(t));
            L1=[t;intx];
            
            if length(dep_prof(:,1))==1;
                L2=[t;dep_prof];
            else
                L2=[t;dep_prof'];
            end
            P = InterX(L1,L2);
            %% if there a multiple, take the 1st of P to represent the point furthest from shore and round to a whole number to index
            if ~isempty(P);
                P=round(P(1,1));
            else
                %If it is empty, that means there are no crosses near 15m
                %and just want to grab the deepest part of the transect
                [~,P]=max(dep_prof);
               
            end   
           I=P;
        end

    else
        %Here the edge is below water, so now if there are any regions
        %above water, limit profile to the furthest edge
        idx=Above(end);
        new_prof=dep_prof;
        new_prof(1:idx-1)=NaN;
        sub_prof=new_prof-h1;
        [Nah ind]=min(abs(sub_prof));
        cut=new_prof(ind);  
        
        
        new_prof(new_prof>cut)=NaN;
        %Now find deepest depth in the remaining profile, which should be
        %closest to 15 in this region
        [Num I]=max(new_prof); %I is the index of the point of interest

        
        if Num<5
           %Here figure out how many blocks there are
            Above=find(dep_prof<0);
            temp=diff(Above);
            temp=find(temp>dif_tol);
            temp=temp+1;
            num_blocks=numel(temp)+1;
           
    
            if num_blocks==1;
                 new_prof=dep_prof;
                 sub_prof=new_prof-h1;
                 [Nah ind]=min(abs(sub_prof));
                 I=ind;

            else
                %Now if there are more than 1 block out there it gets more
                %complicated because need to account for every eventuality
            
                 new_prof=dep_prof;
                 sub_prof=new_prof-h1;
                 [Nah ind]=min(abs(sub_prof));
                 I=ind;

            end
        end
                    
                    
                    
                    
           
    end
else %If the end is above sea level but only that end
    Above=find(dep_prof<0);
    temp=diff(Above);
    temp=find(temp>dif_tol);
    num_blocks=numel(temp)+1;
    
    if num_blocks==1; %essentially just the end above sea level
        if dep_prof(1)>h1 %if its deeper, likely mean on open coast
        sub_prof=dep_prof-h1;
        [Nah ind]=min(abs(sub_prof));
        cut=dep_prof(ind);    
        
        
        dep=find(dep_prof<=cut);
        I=dep(1);

        
        else % if not, likely in a bay and want the closest to 15m as close to shore as possible
            %first see if there is an intersection between the two 
            t=1:length(dep_prof);
            intx=repmat(h1,1,length(t));
            L1=[t;intx];
            if length(dep_prof(:,1))==1;
                L2=[t;dep_prof];
            else
                L2=[t;dep_prof'];
            end
            P = InterX(L1,L2);
            %% if there a multiple, take the first of P to represent the point farthest offshore and round to a whole number to index
            if ~isempty(P);
                P=round(P(1,1));
            else
                %If it is empty, that means there are no crosses near 15m
                %and just want to grab the deepest part of the transect
                [~,P]=max(dep_prof);
               
            end   
           I=P;
        end
    else %There are offshore regions above sea level, want to remove them from consideration if they are less than a certain width
        new_prof=dep_prof;
        temp=diff(Above);
        temp2=find(temp>dif_tol);
        temp2=temp2(end); %only select the point closest to shore
        idx=Above(temp2)+1; %gets index of the wet poiont shoreward of the outcrop
        new_prof(1:idx-1)=NaN;
        
        
        sub_prof=new_prof-h1;
        [Nah ind]=min(abs(sub_prof));
        cut=new_prof(ind); 
        
        new_prof(new_prof>cut)=NaN;
        %Now find deepest depth in the remaining profile, which should be
        %closest to 15 in this region
        [Num I]=max(new_prof); %I is the index of the point of interest
        
        
        
        if Num<5
           %Here figure out how many blocks there are
            Above=find(dep_prof<0);
            temp=diff(Above);
          
            temp=find(temp>dif_tol);
            
            temp=temp+1;
            temp2=temp;
            num_blocks=numel(temp)+1;
            parse_above=1:length(Above);
            for ii=1:num_blocks
                if ii==num_blocks
                    parse_above(temp2(ii-1):end)=ii;
                else
                parse_above(1:temp2(ii))=ii;
                end
            end
           
            if num_blocks==2; %just at the shoreline and offshore
                new_prof=dep_prof;

                new_prof(Above(parse_above==1):end)=NaN;
                sub_prof=new_prof-h1;
                [Nah ind]=min(abs(sub_prof));
                cut=new_prof(ind); 
                if new_prof(1)<h1;
                    I=ind;
                else
                dep=find(new_prof<=cut);
                qz=dep_prof(dep);
                [wz wi]=findpeaks(qz);
                
                if ~isempty(wz);
                [~,ee]=max(wz);
                wz(wi<ee)=[];
                wi(wi<ee)=[];
                
               if numel(wi)==1;
                   I=dep(wi);
               else
                wr=diff(wz);
                
                wa=find(wr>0);
                if ~isempty(wa);
                wa=wa(end);
                wa=wi(wa+1);
               
                I=dep(wa);
                else
                wa=wi(1);
               
                I=dep(wa);
                end
                
                
                
                if dep_prof(I)<10;
                new_prof=dep_prof;
                 sub_prof=new_prof-h1;
                 [Nah ind]=min(abs(sub_prof));
                 I=ind;
                end
                
               end
                else
                  new_prof=dep_prof;

                  new_prof(Above)=NaN;
                  sub_prof=new_prof-h1;
                  [Nah ind]=min(abs(sub_prof));   
                  I=ind;
                end
                
                end
            else
                
                %In this scenario, ignoring the rocks that are in the
                %way of the transect, so just pull
                %the value that is closest to 15m 
                
                %Now if there are more than 1 block out there it gets more
                %complicated because need to account for for every eventuality
                new_prof=dep_prof;
                 sub_prof=new_prof-h1;
                 [Nah ind]=min(abs(sub_prof));
                 I=ind;

            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
end

%So now, I should be the index of the point I want to pull from the
%remaining 



