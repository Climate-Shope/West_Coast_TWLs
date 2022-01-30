%Prep input txt file for use with NOAA's vdatum tool to generate MHW and
%MSL values

NAME_GRD_transects  = 'Douglas'; 

dirout      = ['E:\West_Coast_TWL_Hazards\03_Results\MHW_pre']; if ~exist([dirout],'dir'), mkdir(dirout); end
fname=[NAME_GRD_transects '_MHW.txt'];
fpath=[dirout filesep fname];

%Load transect locations (this will chnage for the different ways the
%transects were calculated)


trans_dirin=['E:\West_Coast_TWL_Hazards\01_Data\Transects\' NAME_GRD_transects '_Transects_v3.mat'];
load(trans_dirin);

eval(['base='  NAME_GRD_transects '_Transects_v3;']);


X=base.X(:,2);
Y=base.Y(:,2);
R=length(X);
Z=zeros(R,1); %<- create a vector of zeros to represent the baseline 
Mat=[X Y Z];

dlmwrite(fpath,Mat,'delimiter','\t','precision',12);




