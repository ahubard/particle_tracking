%% Get xot and  yot for each picture
tic
%% Experimental Info
load('/Users/Aline/Documents/Research/MATLAB/Avalanches/rot_info.mat')
Tr = 2000;
Total_rotation_steps = 1:5:Tr;
Nr = length(Total_rotation_steps);
fnb = '/Users/Aline/Documents/Research/MATLAB/Avalanches/snap_center_1.bin';

%% Ranges of circle parameters
R_range = 617:2:619;
xo_range = 627:1:637;
yo_range = 410:1:420;

%% Saving variables
xot = zeros(1,Nr);
yot = zeros(1,Nr);
Rt = zeros(1,Nr);
peak_H = zeros(1,Nr);
%% main loop
for ii = 1:Nr
    ir = Total_rotation_steps(ii);
    ima = readrot(fnb,info.Nx,info.Ny,ir,0,0); %read image
    ime = edge(ima,'Canny');   %find edges of image
    [iy, ix] = find(ime > 0);
    H = sum(hugh_circle(ime,xo_range,yo_range,R_range),3);
    [mH, I] = max(H(:));
    peak_H(ii) = mH;
    [a, b] = ind2sub(size(H),I);
    aux_xo = xo_range(a);
    aux_yo = yo_range(b);
    ix = ix(iy <300);
    iy = iy(iy<300);
    R = sqrt((iy-aux_yo).^2+(ix-aux_xo).^2);
    ix = ix(R >= (R_range(1)-0.5) & R <= R_range(end)+0.5);
    iy = iy(R >= R_range(1)-0.5 & R <= R_range(end)+0.5);
    circle_parameters = CircleFitByPratt([ix iy]);
    xot(ii) = circle_parameters(1);
    yot(ii) = circle_parameters(2);
    Rt(ii) = circle_parameters(3);
end

%% Find the image points that belonge


%% Fit data to sin, cos

fit_x = fit(Total_rotation_steps',xot','Fourier1');
fit_y = fit(Total_rotation_steps',yot','Fourier1');

toc
% Put use  xot = xo+A_xsin(wth+phi_x) and yot = yo+A_ycos(wth+phi_y)
phi_x = atan(fit_x.a1/fit_x.b1);
phi_y = atan(-fit_y.b1/fit_y.a1);
A_x = fit_x.a1/sin(phi_x);
A_y = fit_y.a1/cos(phi_y);

if(A_x < 0)
    A_x = -A_x;
    phi_x = phi_x + pi;
end

if(A_y < 0)
    A_y = -A_y;
    phi_y = phi_y + pi;
end

file_for_center = 'circle_center_1.mat';
save(file_for_center,'xot','yot','Rt','fit_x','fit_y','phi_x','phi_y','A_x','A_y');
