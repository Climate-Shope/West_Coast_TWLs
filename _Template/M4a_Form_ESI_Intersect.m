%Checks for the intersection of the profile with Envrionmental Sensitivity
%Index and coastal Armoring data for use with runup calculaitons

warning off all 

NAME_GRD_transects  = 'Douglas'; 
addpath('F:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions');


dirout      = ['F:\West_Coast_TWL_Hazards\03_Results\' NAME_GRD_transects '\Profile_ESI']; if ~exist([dirout],'dir'), mkdir(dirout); end

file_esi          = ['F:\West_Coast_TWL_Hazards\01_Data\WashOre_ESI_Shapefiles' filesep 'Wash_Ore_ESI_UTM10_v2.shp'];
file_riprap      = ['F:\West_Coast_TWL_Hazards\01_Data\WashOre_ESI_Shapefiles' filesep 'OR_RipRap_UTM10_v2.shp'];
file_seawall     = ['F:\West_Coast_TWL_Hazards\01_Data\WashOre_ESI_Shapefiles' filesep 'OR_Seawalls_UTM10_v2.shp'];
S=shaperead(file_esi);
T=shaperead(file_riprap);
U=shaperead(file_seawall);
load(['F:\West_Coast_TWL_Hazards\_STEP\_', NAME_GRD_transects,filesep, 'Transects2.mat']);

ESI_Array=[];
g=[];
 for ii=1:length(Transects2.xpf);
     xpf=Transects2.xpf(ii,:);
     ypf=Transects2.ypf(ii,:);
     L1=[xpf;ypf];
     
     if isnan(xpf(1));
           tempx=NaN;
           tempy=NaN;
           esi=NaN;
           Land=NaN;
           Trans=ii;
           
           g={Trans,tempx,tempy,esi,Land,NaN};
           ESI_Array=[ESI_Array;g];
     else
     
     for jj=1:numel(S);
       x=S(jj).X(1:end-1); 
       y=S(jj).Y(1:end-1);
       L2=[x;y];
       
       H1=mean(L1');
       H2=mean(L2');
       
       xcomp=(H1(1)-H2(1))^2;
       ycomp=(H1(2)-H2(2))^2;
       tdistance=sqrt(xcomp+ycomp);
       
       if tdistance > 5000
       continue
       else
       P=InterX(L1,L2);
       
       
       
       if isempty(P)
           continue
       else
           
           
           
           %keep the interseciton that is closest to the
           %offshore point
           xcomp=(P(1,:)-xpf(1)).^2;
           ycomp=(P(2,:)-ypf(1)).^2;
           distance=sqrt(xcomp+ycomp);
           [~,I]=min(distance);
           
           
           
           tempx=P(1,I);
           tempy=P(2,I);
           esi=S(jj).ESI;
           k=strfind(esi,'/');
           
           if ~isempty(k);
           esi(k(1):end)=[];
           end 
           
           Land=S(jj).LANDWARD_S;
           Trans=ii;
           
           g={Trans,tempx,tempy,esi,Land,NaN};
           ESI_Array=[ESI_Array;g];
           
           
       end
       end
     end
     
          if exist('g','var')
         if isempty(g);
         tempx=NaN;
          tempy=NaN;
           esi=NaN;
           Land=NaN;
           Trans=ii;
           
           g={Trans,tempx,tempy,esi,Land,NaN};
           ESI_Array=[ESI_Array;g];
         end
     end
  if ~exist('Trans','var')
         tempx=NaN;
          tempy=NaN;
           esi=NaN;
           Land=NaN;
           Trans=ii;
           
           g={Trans,tempx,tempy,esi,Land,NaN};
           ESI_Array=[ESI_Array;g];
     end
     clear g Trans    
     
     
     
     end
 end
      
 
 
 
 %% Now write over the structures 
 for ii=1:length(Transects2.xpf);
     xpf=Transects2.xpf(ii,:);
     ypf=Transects2.ypf(ii,:);
     L1=[xpf;ypf];
     
     if isnan(xpf(1));
        continue
     else
     
     for jj=1:numel(T);
       x=T(jj).X(1:end-1); 
       y=T(jj).Y(1:end-1);
       L2=[x;y];
       
       H1=mean(L1');
       H2=mean(L2');
       
       xcomp=(H1(1)-H2(1))^2;
       ycomp=(H1(2)-H2(2))^2;
       tdistance=sqrt(xcomp+ycomp);
       
       if tdistance > 5000
       continue
       else
       P=InterX(L1,L2);
       
       
       
       if isempty(P)
           continue
       else
           
           
           
         %Fix for Riprap
           
            esi='6B';
  
            nn=cell2mat(ESI_Array(:,1));
           mm=find(ismember(nn,ii));
           
           for ttt=1:length(mm);
           %ESI_Array{mm(ttt),6}=T(jj).STRUCTURE;
           ESI_Array{mm(ttt),6}='riprap';
           
           Trans=mm(ttt);
           if strcmp(ESI_Array{Trans,4},esi)==0;
           ESI_Array{Trans,4}=esi;
           end
           end
           
           
       end
       end
     end
     

     
     
     
     end
 end
 
 
 %% Now write for the detached seawalls
 %Detached seawalls were identified spearately in GIS and converted to
 %shapefiles
 SeaIntersect=[];
 for ii=1:length(Transects2.xpf);
     xpf=Transects2.xpf(ii,:);
     ypf=Transects2.ypf(ii,:);
     L1=[xpf;ypf];
     
     if isnan(xpf(1));
        continue
     else
     
     for jj=1:numel(U);
       x=U(jj).X(1:end-1); 
       y=U(jj).Y(1:end-1);
       L2=[x;y];
       
       H1=mean(L1');
       H2=mean(L2');
       
       xcomp=(H1(1)-H2(1))^2;
       ycomp=(H1(2)-H2(2))^2;
       tdistance=sqrt(xcomp+ycomp);
       
       if tdistance > 5000
       continue
       else
       P=InterX(L1,L2);
       
       
       
       if isempty(P)
           continue
       else
           
           %keep the intersection that is closest to the
           %offshore point
           xcomp=(P(1,:)-xpf(1)).^2;
           ycomp=(P(2,:)-ypf(1)).^2;
           distance=sqrt(xcomp+ycomp);
           [~,I]=min(distance);
           
           Int.X=P(1,I);
           Int.Y=P(2,I);
           temp=[ii Int.X Int.Y]; 
           SeaIntersect=[SeaIntersect;temp];
           
           
           %Fix for Seawall

           esi='1B';
 
           nn=cell2mat(ESI_Array(:,1));
           mm=find(ismember(nn,ii));
           
           for ttt=1:length(mm);
           ESI_Array{mm(ttt),6}='seawall';
           
           
           Trans=mm(ttt);
           if strcmp(ESI_Array{Trans,4},esi)==0;
           ESI_Array{Trans,4}=esi;
           end
           end
           
           
       end
       end
     end
     

     
     
     
     end
 end
 
 
 
save([dirout,filesep,NAME_GRD_transects,'_ESI_Array'],'ESI_Array');   
save([dirout,filesep,NAME_GRD_transects,'_SeaIntersect'],'SeaIntersect');        
       
       