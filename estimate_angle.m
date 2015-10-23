function [theta, slope, xy_mass, ii_surface, bindex] = estimate_angle(PX,PY,D,yo,Nb_bins,giveCM,NT)

%Finds angle of surface of the pile. giveCM = 0 returns only the angle
% giveCM = 1; Finds the center of mass on each bin in the surface for t -1;



binsize = D;

Ly =400;
if(~exist('NT','var'))
    NT = size(PX,2);
end
theta = zeros(1,NT); 
slope = zeros(2,NT);



        

for t = 1:NT
    
    
    x = PX(:,t);
    y = PY(:,t);
    
    %Find y boundaries for surface of the pile
    y_top = min(y);
    y_bottom = yo + max(4*D,yo-y_top);
    
    ii_surface = find(y < y_bottom);
    x = x(ii_surface);
    y = y(ii_surface);
    
    % put data in bins
    [xbin, ix] = sort(ceil(x/binsize));
    bindex = [0 ;find(diff(xbin)>0)];
    y = y(ix);
    x = x(ix);
    
    % Angle variable
    yangle = zeros(1,length(bindex)-1);
    xangle = yangle;
    
    if((giveCM > 0) && (t ==1))
        ii_surface = ii_surface(ix);
        %Center of mass bin
        xy_mass = zeros(2,Nb_bins);
    end
    
    for bin = 1:length(bindex)-1;
        yangle(bin) = min(y(bindex(bin)+1:bindex(bin+1)));
        xangle(bin) = binsize*bin;
        if ((giveCM > 0) && (t ==1))
            total_m = length(bindex(bin)+1:bindex(bin+1));
            xy_mass(1,bin) = sum(x(bindex(bin)+1:bindex(bin+1)))/total_m;
            xy_mass(2,bin) = sum(y(bindex(bin)+1:bindex(bin+1)))/total_m;    
        end
    end
    
    top_layer = find(yangle>0);
    xangle = xangle(top_layer);
    yangle = Ly - yangle(top_layer);
    
    slope(:,t) = polyfit(xangle,yangle,1);
    theta(t) = atan(slope(1,t))*180/pi;
end
    
    


