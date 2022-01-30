%After Vdatum calculation, extracts MSL levels from Vdatum output txt file

NAME_GRD_transects  = 'Douglas'; 

dirout      = ['E:\West_Coast_TWL_Hazards\03_Results\MSL']; if ~exist([dirout],'dir'), mkdir(dirout); end
fname=[NAME_GRD_transects '_MSL.mat'];
fpath=[dirout filesep fname];

%Load transect locations (this will change for the different ways the
%transects were calculated)


dirin=['E:\West_Coast_TWL_Hazards\03_Results\MSL_post'];
fnamein=[NAME_GRD_transects '_MHW.txt'];

fpath2=[dirin filesep fnamein];

M=dlmread(fpath2,'\t',1,0);

MSL=M(:,6)*-1;%the datum is negative relative to the 0 value

MSL(MSL>10)=NaN;
num=[1:length(MSL)]';
numblank=num;
numblank(isnan(MSL))=[];
MSL(isnan(MSL))=[];
%Interp if Vdatum returned NaN values
Vout = interp1(numblank,MSL,num);
MSL=Vout;


X=M(:,4);
Y=M(:,5);

MSL=[X Y MSL];

save(fpath,'MSL');





