function [independentData,dailyMaxValues,threshstart] = createIndependentEventgeneric(data,serialTime,valmin,delta_t)
%  CREATEINDEPENDENTEVENTTS creates a time series of independent events
%
% [independentWaveNTR,dailyMaxValues] =
% createIndependentEventTS(inpath,filename,valmin,delta_t,outpath,station,dateid)
% INPUT 
% inpath:  folder containing combinedWaveTide.mat buoy/tide gauge input 
% files created from function createCombinedWaveTideTS.m
% 
%             example - 'D:/matlab/generic/inpath/'         
% filename: Name of file in above inpath, should be CombinedWaveTide.mat
% outpath:  folder where progam will save Matlab output
%               example - 'D:/matlab/generic/outpath'
% stationid: numeric; station name assigned by user 
% dateid: string; date this file was created (yyyymmdd)
% OPTIONAL
% valmin: The minimum value of nonNan wave heights in a day (usually 3)
% delta_t: timestep to make sure events are independent (72 hours or 3
% days)
% OUTPUT
% 
% (1) Creates variable 'independentData' and saves to outpath.
%       The size of this variable is Nx2, where N is the number of
%       daily measurements.
%           Column 1: serial time
%           Column 2: variable
% (2) Creates variable 'threshstart' which is the threshold to use for the GPD 


% HISTORY
% 2012: Developed by Katy Serafin(OSU)
% Last updated by K.Serafin 1/3/2013

%% If optional values are empty, 
% check for a save id, if it doesn't exist, save
if ~exist('valmin'), valmin = 3; end
if ~exist('delta_t'), delta_t = 3; end


date_vector= datevec(serialTime);

year = date_vector(:,1);
month = date_vector(:,2);
day = date_vector(:,3);
hour = date_vector(:,4);
minute = date_vector(:,5);


% Remove any data that starts earlier than wave height measurements begin
clop=find(isnan(data)==0);
ind=clop(1);
out = 1:ind-1;
out=out';

year(out) = [];
month(out) = [];
day(out) = [];
hour(out) = [];
minute(out) = [];
data(out) = [];
serialTime(out) = [];
%% Generate a timeseries of daily maximum values; need a minimum of valmin
%% values in a day

len = floor(serialTime(end))-floor(serialTime(1));      % number of days
waveday = floor(serialTime(1)):1:floor(serialTime(end));% day 1 is our first
                                                        % day of measurements
cnt1 = 0;
for ii = 1:length(waveday);
    ind1 = find(floor(serialTime) == waveday(ii));
    if length(find(data(ind1)>0)) > valmin;
        [day_max(ii),I] = max(data(ind1)); % find max wave height for each waveday
        ind_max(ii) = ind1(I);
        time_max(ii) = waveday(ii);
   elseif isempty(ind1)
        cnt1 = cnt1+1;
        day_max(ii) = NaN;
        ind_max(ii) = NaN;
        other_ind(cnt1) = NaN;
        time_max(ii) = waveday(ii);
    else
        cnt1 = cnt1+1;
        day_max(ii) = NaN;
        ind_max(ii) = ind1(1);
        time_max(ii) = waveday(ii);
        other_ind(cnt1) = ii;
    end
    clear ind1
end

%% Create a daily maximum time series
day_fin = serialTime(ind_max(find(isnan(ind_max)<1)));
% Get rid of dates where there weren't enough data points (valmin) to be
% valid
DATA = data;
%DATA(ind_max(other_ind)) = NaN;
height_fin = DATA(ind_max(find(isnan(ind_max)<1)));

% Create a structure of these to save for later
dailyMaxValues.day_fin = day_fin;
dailyMaxValues.height_fin = height_fin;

% delta t filter; gets rid of data +- delta t away to eliminate any
% dependency

% Set up a vector for all parameters
independentWaveNTR = nan(length(day_fin),2);
independentWaveNTR(:,1) = day_fin;


ind = find(day_fin);
[hkeep,I] = max(height_fin(ind));
dkeep = day_fin(I);

independentWaveNTR(I,2) = height_fin(I);



% Now isolate all events and only include those that happen delta_t days
% apart
while length(ind) > delta_t
    ind2 = ind(find(abs(day_fin(ind)-day_fin(ind(I)))>=delta_t));
    [hdum,I] = max(height_fin(ind2));
    hkeep = [hkeep;hdum];
    dkeep = [dkeep;day_fin(ind2(I))];
    independentWaveNTR(ind2(I),2) = height_fin(ind2(I));
    ind = ind2;
end

clear ind

%% Save starting estimates for threshold as well....
[f,x] = ecdf(data);
threshstart.data=interp1(f,x,0.95);

% save for later use
independentData = independentWaveNTR;

end

