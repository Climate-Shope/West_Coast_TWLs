%Code to smooth the shoreline of the complex Norhwest Shorleines so there
%are not as many errant transects (still may require editing by hand)

clear all; run('E:\West_Coast_TWL_Hazards\_STEP\_STEP_Functions\init_config.m') 

NAME_GRD_transects  = 'Douglas'; 

%Provide a user defined smoothing window
window=2;
window_name=10*(2*(window))+10;


%Load relevant directories and shoreline position information
run('run_init_directories.m')
load([NAME_GRD_transects '_Shoreline.mat']);
eval(['Points.X=' NAME_GRD_transects '_Shoreline.X;']);
eval(['Points.Y=' NAME_GRD_transects '_Shoreline.Y;']);
Points=[Points.X', Points.Y'];
Points=unique(Points,'rows','stable');


%Step1: get total shoreline distance
Xcomp=diff(Points(:,1)).^2;
Ycomp=diff(Points(:,2)).^2;
Distance=sqrt(Xcomp+Ycomp);

%Calculate the cumulative distance vector
CumDistance=cumsum(Distance);
CumDistance=[0;CumDistance];

numpts=round(CumDistance(end)/10);
vert=linspace(0,CumDistance(end),numpts);

test=[CumDistance,Points(:,1),Points(:,2)];
test=unique(test,'rows','legacy');
X=interp1(test(:,1),test(:,2),vert);
Y=interp1(test(:,1),test(:,3),vert);

X=runmean(X,window); %<-- Results in a 30m moving average window
Y=runmean(Y,window); %<-- Results in a 30m moving average window


%Step 2: Re-interpolate the distances and recast points every 10m
Xcomp=diff(X).^2;
Ycomp=diff(Y).^2;
Distance=sqrt(Xcomp+Ycomp);
%Calculate the cumulative distance vector
CumDistance=cumsum(Distance);
CumDistance=[0,CumDistance];
numpts=round(CumDistance(end)/10);
vert=linspace(0,CumDistance(end),numpts);
test=[CumDistance',X',Y'];
test=unique(test,'rows','legacy');
X=interp1(test(:,1),test(:,2),vert);
Y=interp1(test(:,1),test(:,3),vert);


%Step 3: Now save the outputs as a new shoreline

Shoreline.X=X;
Shoreline.Y=Y;

name=[NAME_GRD_transects '_Shoreline_' num2str(window_name) 'm_Smooth'];

eval([name '=Shoreline;']);

save(name,name);








