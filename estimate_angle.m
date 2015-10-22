function [theta, slope]= estimate_angle(PX,PY,D, NT,R,xo,yo)

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
    y_bottom = yo + max(2.5*D,yo-y_top);
    
    ii_surface = find(y < y_bottom);
    x = x(ii_surface);
    y = y(ii_surface);
    
    % put data in bins
    [x, ix] = sort(ceil(x/binsize));
    bindex = [0 ;find(diff(x)>0)];
    y = y(ix);
    ii_surface = ii_surface(ix);
    
    % Angle variable
    yangle = zeros(1,length(bindex)-1);
    xangle = yangle;
    
    %Center of mass bin
    x_mass = zeros(1,length(bindex)-1);
    y_mass = zeros(1,length(bindex)-1);
    
    for bin = 1:length(bindex)-1;
        yangle(bin) = min(y(bindex(bin)+1:bindex(bin+1)));
        xangle(bin) = binsize*bin;
        total_m = length(bindex(bin)+1:bindex(bin+1));
        x_mass(bin) = sum(x(bindex(bin)+1:bindex(bin+1)))/total_m;
        y_mass(bin) = sum(y(bindex(bin)+1:bindex(bin+1)))/total_m;    
    end
    
    top_layer = find(yangle>0);
    xangle = xangle(top_layer);
    yangle = Ly - yangle(top_layer);
    
    slope(:,t) = polyfit(xangle,yangle,1);
    theta(t) = atan(slope(1,t))*180/pi;
end
    
    


