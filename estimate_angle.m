function [theta, slope]= estimate_angle(PX,PY)

binsize = 10;
Ly =400;
NT = size(PX,2);
theta = zeros(1,NT); 
slope = zeros(2,NT);

for t = 1:NT
    x = PX(:,t);
    y = PY(:,t);
    
    [x, ix] = sort(ceil(x/binsize));
    bindex = [0 ;find(diff(x)>0)];
    y = y(ix);
    ya = zeros(1,length(bindex)-1);
    xa = ya;
    
    for bin = 1:length(bindex)-1;
        ya(bin) = min(y(bindex(bin)+1:bindex(bin+1)));
        xa(bin) = binsize*bin;
    end
    
    xa = xa(ya>0);
    ya = Ly - ya(ya>0);
    
    slope(:,t) = polyfit(xa,ya,1);
    theta(t) = atan(slope(1,t))*180/pi;
end
    
    


