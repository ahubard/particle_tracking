function [slope, theta]= estimate_angle(PX,PY)

binsize = 10;
Ly =400;



    
    
    
    x = PX(PX>0);
    y = PY(PX>0);
    [x ix] = sort(ceil(x/binsize));
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
    
    slope = polyfit(xa,ya,1);
    theta = atan(slope(1))*180/pi;
    
    


