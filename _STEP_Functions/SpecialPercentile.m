function [Y] = SpecialPercentile(arr, pct)
%Generate all of the whole values by fitting
%a curve to the arr and then interpolating to whichever point is there using
%the polyfit function

arr(isnan(arr))=[];

if pct==100;
    Y=max(arr);
else

pct=pct/100;
num=1:length(arr);
num=num-1;
arr=sort(arr)';
pct=num(end)*pct; %represents the points in the arr thats closest

%find the closest 2 points
temp=abs(pct-num);
[~,I]=min(temp);
one=I;
temp(one)=inf;
[~,I]=min(temp);
two=I;
pct=pct+1; %To adjust the percent to an index between 1 and length(arr) 
%now organize in ascending order;
points=sort([one,two]);

%interpolate between the two
Y=interp1(points,arr(points),pct);
end

end