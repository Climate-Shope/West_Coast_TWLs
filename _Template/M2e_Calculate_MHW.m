%After Vdatum calculation, extracts MHW levels from Vdatum output txt file


NAME_GRD_transects  = 'Douglas'; 

dirout      = ['E:\West_Coast_TWL_Hazards\03_Results\MHW']; if ~exist([dirout],'dir'), mkdir(dirout); end
fname=[NAME_GRD_transects '_MHW.mat'];
fpath=[dirout filesep fname];

%Load transect locations (this will change for the different ways the
%transects were calculated)


dirin=['E:\West_Coast_TWL_Hazards\03_Results\MHW_post'];
fnamein=[NAME_GRD_transects '_MHW.txt'];

fpath2=[dirin filesep fnamein];

M=dlmread(fpath2,'\t',1,0);

MHW=M(:,6)*-1; %the datum is negative relative to the 0 value

MHW(MHW>10)=NaN;
num=[1:length(MHW)]';
numblank=num;
numblank(isnan(MHW))=[];
MHW(isnan(MHW))=[];
%Interp if Vdatum returned NaN values
Vout = interp1(numblank,MHW,num);
MHW=Vout;

X=M(:,4);
Y=M(:,5);

MHW=[X Y MHW];

save(fpath,'MHW');





