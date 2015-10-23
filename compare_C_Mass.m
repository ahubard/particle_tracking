function [ d_mass] = compare_C_Mass(x,y,bindex,xy_mass)

%Finds angle of surface of the pile. giveCM = 0 returns only the angle
% giveCM = 1; Finds the center of mass on each bin in the surface for t -1;








 d_mass = zeros(2,length(bindex)-1);
    for bin = 1:length(bindex)-1;
      
            total_m = length(bindex(bin)+1:bindex(bin+1));
            d_mass(1,bin) = sum(x(bindex(bin)+1:bindex(bin+1)))/total_m;
            d_mass(2,bin) = sum(y(bindex(bin)+1:bindex(bin+1)))/total_m;    
        
    end
    
   d_mass = sqrt(sum((xy_mass - d_mass).^2));