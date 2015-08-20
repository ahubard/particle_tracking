function [imeanangle fmeanangle]= estimate_angle(En)

datafile = sprintf('/aline/rotdrum/o%i/AVALANCHE_%i.mat',En,En);
load(datafile,'avaini','avaendt','avainit');
binsize = 10;
Numangles = length(avaini);
initialangle = zeros(1,Numangles)-1;
finalangle = zeros(1,Numangles)-1;
Ly =400;


for nb = 1:Numangles
    
    filenb = avaini(nb);
    tini = avainit(nb);
    tend = avaendt(nb);
    
    posfile = sprintf('/aline/rotdrum/o%i/Tracked_%i%03i.mat',En,En,filenb);
    load(posfile,'PX','PY');
    
    
    
    x = PX(PX(:,tini)>0,tini);
    y = PY(PX(:,tini)>0,tini);
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
    initialangle(nb) = atan(slope(1))*180/pi;
    
    
    
    x = PX(PX(:,tend)>0,tend);
    y = PY(PX(:,tend)>0,tend);
    
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
    finalangle(nb) = atan(slope(1))*180/pi;
end

imeanangle = mean(initialangle+29);
fmeanangle = mean(finalangle+29);
saveavafile = sprintf('/aline/rotdrum/o%i/AVALANCHE_%i.mat',En,En);
    
save(saveavafile,'initialangle','finalangle','-append'); 