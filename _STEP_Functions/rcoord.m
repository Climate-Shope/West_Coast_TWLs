function [u,v] = rcoord(r, az)
% Polar to rectangular coordinate conversion, geographic convention
%
% function [u,v] = xycoord(speed, sdir)
% function [x,y] = xycoord(r, az)
%
% Converts from speed, azimuth polar coordinates (0-360 degrees)
% to rectangular coordinates. (Note: use of dir for direction conflicts
% with DOS/MATLAB dir command!)

% Jessie Lacy, USGS
% May 26, 2010
% last modified July 20, 2010

theta = (90-az)*pi/180;     %radians, trig convention

u = r .* cos(theta);
v = r .* sin(theta);
return

