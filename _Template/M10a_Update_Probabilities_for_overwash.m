% Script that fixes small probability issues in the Final probabilities for 
% pCOI

%First describe the appropriate directory 
Dirout='E:/West_Coast_TWL_Hazards/03_Results/Return_Levels_updated_03_17_2021';
if ~exist(Dirout,'dir'), mkdir(Dirout); end 
Dirin='E:/West_Coast_TWL_Hazards/03_Results/Return_Levels_updated_03_16_2020';
cd(Dirout);
% Return_Periods
periods={'1','2','5','10','20','25','50','100','250','500'};

%Load in the txt file 

for ii=1:length(periods)
  %Fixes the probabilites where there are negative overwash probabilities
  a=['RP' periods{ii} '_Westcoast.mat']
  b=load([Dirin filesep a]);
  eval(['b=b.RP' periods{ii} '_Westcoast;']);
  
  c=find(b(1:end,22)<0);
  b(c,21)=b(c,21)+b(c,22);
  b(c,22)=0;
  
  c=find(b(1:end,26)<0);
  b(c,25)=b(c,25)+b(c,26);
  b(c,26)=0;
  
  
  c=find(b(1:end,21)<0 & b(1:end,21)>-0.0001);
  b(c,21)=0;
 
  
  
  %Next need to do a text write with headers in the dir out field 
  filename=[Dirout filesep 'RP' periods{ii} '_Westcoast.txt']
  dlmwrite (filename, b, ",")
  
  clear a b c
end
