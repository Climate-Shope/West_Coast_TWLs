% Script that does an outlier analysis on the alongshore TWL values and
% removes errant profiles

%First describe the appropriate directory 
Dirout='E:/West_Coast_TWL_Hazards/03_Results/Return_Levels_updated_04_06_2021';
if ~exist(Dirout,'dir'), mkdir(Dirout); end 
Dirin='E:/West_Coast_TWL_Hazards/03_Results/Return_Levels_updated_03_17_2021';
cd(Dirout);
% Retun_Periods
periods={'1','2','5','10','20','25','50','100','250','500'};

%Load in the txt file 

for ii=1:length(periods)
  %Fixes the probabilites where there are negative overwash probabilities
  a=[Dirin filesep 'RP' periods{ii} '_Westcoast.txt'];
  b=dlmread(a);
  cc=[];
  for jj=1:length(b(:,1));
    %if the TWL is < MHW, need to remove point
    if b(jj,5)< b(jj,17)
      cc=[cc,jj];
    else
    
    
    %isolate the point in question
    aa=[b(jj,1),b(jj,2)];
    aa_twl=b(jj,5);
    
    %calculate the dist between this point and all of the rest
    xx=b(:,2);
    yy=b(:,1);
    tt=b(:,5);
    
    
    dx=(aa(2)-xx).^2;
    dy=(aa(1)-yy).^2;
    d=(dx+dy).^(1/2);
    
    dind=find(d<0.01);
    
    
    %now just use the new vector to get the 
    tt=tt(dind);
    d=d(dind);
    
    dind=find(d~=0);
    tt=tt(dind);
    
    %now get the mean and standard deviation of the range
    t_avg=mean(tt);
    t_std=std(tt);
    t_std3=t_std*3;
    cutoff=t_avg+t_std3;
    
    %now log if the point needs to be removed, if it is more than 3 std of
    %the regional mean 
    
    if aa_twl>cutoff
      cc=[cc,jj];
    end
    
    clear aa aa_twl xx yy tt dx dy d dind t_avg t_std t_std3 cutoff 
    
    endif
   end
    
 %now use the cc matrix to remove the rows with issues and write new file 
 b(cc,5)=NaN; 
 b(cc,6)=NaN; 
 b(cc,7)=NaN; 
 b(cc,8)=NaN; 
 b(cc,9)=NaN; 
 b(cc,10)=NaN; 
 b(cc,20)=NaN; 
 b(cc,21)=NaN; 
 b(cc,22)=NaN;
 b(cc,23)=NaN;
 b(cc,24)=NaN;
 b(cc,25)=NaN;
 b(cc,26)=NaN;
 b(cc,27)=NaN;



 filename=[Dirout filesep 'RP' periods{ii} '_Westcoast.txt'];
 dlmwrite (filename, b, ",")   
 clear b cc
   
    
    
  end
  
  


  
end
