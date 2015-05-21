function mask = circlesmask(px,py,NX,NY,Np,R)
% Creates  an image. Image size is NX*NY with circles of Radius R 
% centered ins px,py. Value inside circle is zero, outside one. 

mask = ones(NY,NX);

[localgridy  localgridx] = ndgrid(-R:R,-R:R);
circle = abs(localgridx+1i*localgridy)>R;

xmin = min(px);
xmax = max(px);
dx = xmin-R-1;
Lx = length(xmin:xmax+R*2);

ymin = min(py);
ymax =max(py);
dy =ymin-R-1;
Ly =length(ymin:ymax+R*2);

submask = ones(Ly,Lx);

for np = 1:Np
    submask((-R:R)+py(np)-dy,(-R:R)+px(np)-dx) = submask((-R:R)+py(np)-dy,(-R:R)+px(np)-dx).*circle;
end


mask(max(dy,0)+1:min(Ly+dy,NY),max(dx,0)+1:min(Lx+dx,NX)) = ...
    submask(max(dy,0)+1-dy:min(Ly+dy,NY)-dy,max(dx,0)+1-dx:min(Lx+dx,NX)-dx);