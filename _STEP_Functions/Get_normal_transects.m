function [TransectX, TransectY] = Get_normal_transects(xi,yi,r_on,r_off,num_points);

% Function to obtain shoreline angle for each point in a shoreline vector 
% and then calculate the shore normal transect at each point 

%xi and yi are the coordinates of regularly spaced alongshore points
%r_on is the desired onshore extent of the transect in m
%r_off is the desired offshore extent of the transect in m
%num_points are the points on either side fo the transect used to cal is
%angle. 1 would be one point on either side. 2 would be based on the vector
%between the second point on either side. 

% xi and yi vectors should represent points progresisng south to north 

% TransectX and TransectY have the format of having coordinates that
% represent [offshore, shoreline, onshore] points along the transect

Len=length(xi);

TransectX=[];
TransectY=[];

%%% Start with Point 1
ii=1;
U=xi(ii+num_points)-xi(ii); % Only using the point itself and next in series to get angle because it is first
V=yi(ii+num_points)-yi(ii);
[r,az] = pcoord(U,V);
az_hi= az+90;
az_low=az-90;
[u_hi,v_hi] = rcoord(r_on, az_hi); % Onshore
[u_low,v_low] = rcoord(r_off, az_low); %Offshore
Transectxi = [(xi(ii)+u_low) xi(ii) (xi(ii)+u_hi)];
Transectyi = [(yi(ii)+v_low) yi(ii) (yi(ii)+v_hi)];
TransectX=[TransectX;Transectxi];
TransectY=[TransectY;Transectyi];

%Loop through 2:n-1 points using 2 points on either side to calc shoreline
%angle
for ii=2:(Len-1);
    
     U=xi(ii+num_points)-xi(ii-num_points);
     V=yi(ii+num_points)-yi(ii-num_points);

    
    [r,az] = pcoord(U,V);
    
    az_hi= az+90;
    az_low=az-90;
    
    
    [u_hi,v_hi] = rcoord(r_on, az_hi); % Onshore
    [u_low,v_low] = rcoord(r_off, az_low); %Offshore
    
    % from here need to add these changes to the intial point

    Transectxi = [(xi(ii)+u_low) xi(ii) (xi(ii)+u_hi)];
    Transectyi = [(yi(ii)+v_low) yi(ii) (yi(ii)+v_hi)];
    
    TransectX=[TransectX;Transectxi];
    TransectY=[TransectY;Transectyi];
    
    
end

% Finish with endpoint that will ust use the preceding point to calc angle

ii=Len;
U=xi(ii)-xi(ii-1); % Only using the point itself and prior in series to get angle because it is last
V=yi(ii)-yi(ii-1);
[r,az] = pcoord(U,V);
az_hi= az+90;
az_low=az-90;
[u_hi,v_hi] = rcoord(r_on, az_hi); % Onshore
[u_low,v_low] = rcoord(r_off, az_low); %Offshore
Transectxi = [(xi(ii)+u_low) xi(ii) (xi(ii)+u_hi)];
Transectyi = [(yi(ii)+v_low) yi(ii) (yi(ii)+v_hi)];
TransectX=[TransectX;Transectxi];
TransectY=[TransectY;Transectyi];

%%% Form here, Have transect coordinates that can be used to generate
%%% shapefile. 





    


