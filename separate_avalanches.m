function  Nb_over_boundary = separate_avalanches(filedirectory, En, D, PX,PY)
%% Gives the number of particle that hit the boundary in an avalanche. 
% to separate avalanches afected by the boundaries. 
%% Get boundary info
file_for_boundary = sprintf('%sAvanonestep%i.mat',filedirectory,En);
%boundary is given by circle of radius R centered at (xo,yo),mk gives size
%of the image.
load(file_for_boundary,'mk','R','xo','yo'); 
[Ly, Lx] = size(mk); 
NT = size(PX,2);
%% Left boundary
ydiff2_max = max((1-yo)^2,(Ly-yo)^2);
x_boundary_cuttof = xo-sqrt(R^2-ydiff2_max)+D/2;
Nb_over_boundary = [0 0];
%% Find drop rim point.
% displacement_file = sprintf('Displacement_%i.mat',filenb);
% load(displacement_file,'PX','PY');
ii_near_boundary = find(PX(:,1)<x_boundary_cuttof);
r_near_boundary = (PX(ii_near_boundary,1)-xo).^2 + (PY(ii_near_boundary,1)-yo).^2;
ii_near_boundary = ii_near_boundary(r_near_boundary > (R-.75*D)^2);

if (isempty(ii_near_boundary))
    Nb_over_boundary(1) = 0;
else                                                   
    
    [y_top_layer,  ii_top_layer] = min(PY(ii_near_boundary,1));
    x_top_layer = PX(ii_near_boundary(ii_top_layer),1);
    slope_radial = (y_top_layer-yo)/(x_top_layer-xo);
    slope_boundary = (slope_radial+D/R)/(1-slope_radial*D/R); %tan(D/R)=D/R
    x_boundary = xo-sqrt(R^2/(1+slope_boundary^2));
    y_boundary = slope_boundary*(x_boundary-xo)+yo;
    
    %% Check for particles touching this point in subsequenbt avalanche.
    ii_near_boundary = find(PY(:) < y_boundary);
    ii_near_boundary = ii_near_boundary(PX(ii_near_boundary) < x_boundary+D);
    go_over_boundary = ((R-D/3)^2- ((PX(ii_near_boundary)-xo).^2 + ...
        (PY(ii_near_boundary)-yo).^2) < 0);
    Nb_over_boundary(1) = sum(go_over_boundary);
end

%% Right boundary
ydiff2_max = max((1-yo)^2,(Ly-yo)^2);
x_boundary_cuttof = xo+sqrt(R^2-ydiff2_max)-D/2;


%% Find top point
% displacement_file = sprintf('Displacement_%i.mat',filenb);
% load(displacement_file,'PX','PY');
ii_near_boundary = find(PX(:,1) > x_boundary_cuttof);
r_near_boundary = (PX(ii_near_boundary,1)-xo).^2 + (PY(ii_near_boundary,1)-yo).^2;
ii_near_boundary = ii_near_boundary(r_near_boundary < (R-.75*D)^2);

if (isempty(ii_near_boundary))
    Nb_over_boundary(2) = Nb_over_boundary + 0;
else                                                   
    
    [y_top_layer,  ii_top_layer] = min(PY(ii_near_boundary,1));
    x_top_layer = PX(ii_near_boundary(ii_top_layer),1);
    slope_radial = (y_top_layer-yo)/(x_top_layer-xo);
    slope_boundary = (slope_radial+D/R)/(1-slope_radial*D/R); %tan(D/R)=D/R
    x_boundary = -xo+sqrt(R^2/(1+slope_boundary^2));
    y_boundary = slope_boundary*(x_boundary-xo)+yo;
   
    %% Check for particles close to this point in the first frame that move after. 
    ii_near_boundary = find(abs((PY(:,1) - y_boundary)<1.5*D));
    ii_near_boundary = ii_near_boundary(PX(ii_near_boundary,1) > (x_boundary-D));
    dr_total_near_boundary = (PX(ii_near_boundary,1)-PX(ii_near_boundary,NT)).^2 + ...
        (PY(ii_near_boundary,1)-PY(ii_near_boundary,NT)).^2;
    go_over_boundary = (dr_total_near_boundary > D);
    Nb_over_boundary(2) = Nb_over_boundary(2) + sum(go_over_boundary);
end


    










