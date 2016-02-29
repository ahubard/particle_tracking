function   H = hugh_circle(b_ima, xo_range,yo_range,R_range)

%Finds the hugh transform for a binary image b_ima for circles of radii in
%the R-range with center in the xo_range yo_range with rad.

% Input xo_range, yo_range and R_range are vectors of positive numbers. 
% Output H is a 3dimensional object of the size given by the input vectors
% where 

lx = length(xo_range);
ly = length(yo_range);
lR = length(R_range);

dR = (R_range(lR) - R_range(1))/lR;

H = zeros(lx,ly,lR);
[by, bx] = find(b_ima);
bx = bx(by < 500);
by = by(by < 500);
for ix = 1:lx
    xo = xo_range(ix);
    dx2 = (bx - xo).^2;
    for iy = 1:ly
        yo = yo_range(iy);
        dy2 = (by - yo).^2;
        R = sqrt(dx2 + dy2);
        r = ceil((R - R_range(1))/dR);
        for ir = 1:lR
            H(ix,iy,ir) = H(ix,iy,ir) + sum(r == ir);
        end
    end
end

    
        
        
        
