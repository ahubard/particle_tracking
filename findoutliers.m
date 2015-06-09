function [pxs pys Npf]= findoutliers(folder,En,n)

fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,n);
load(fn,'pxs','pys','Npf');
binsize = 20;
nframes = length(Npf);
D = 10;
nparticles = size(pxs,1);

for nb = 1:nframes
    
    x = pxs(1:Npf(nb),nb);
    y = pys(1:Npf(nb),nb);
    [xs ix] = sort(ceil(x/binsize));
    bindex = [0 ;find(diff(xs)>0)];
    ys = y(ix);
    ya = zeros(1,length(bindex)-1);
    xa = ya;
    for bin = 1:length(bindex)-1;
        ya(bin) = min(ys(bindex(bin)+1:bindex(bin+1)));
        xa(bin) = binsize*(bin + xs(1)-1);
    end
    
   linearfit = polyfit(xa,ya,1);
   nvector = [-linearfit(1) 1]/(sqrt(linearfit(1)^2+1)); %normalvector to bestfit line
   lo = linearfit(2)*nvector(2); %Distance of bestfitline to origin. 
   dpointtoline = x*nvector(1) + y*nvector(2) - lo;
   insiders = (dpointtoline > -D*3); %Points that are in the region of interest
   
   %Get rid of outliers
   Npf(nb) = sum(insiders);
   pxs(1:Npf(nb),nb) = pxs(insiders,nb);
   pxs(Npf(nb)+1:nparticles,nb) = 0;
   pys(1:Npf(nb),nb) = pys(insiders,nb);
   pys(Npf(nb)+1:nparticles,nb) = 0;
   
   
end

delete(fn);
    
save(fn,'pxs','pys','Npf'); 