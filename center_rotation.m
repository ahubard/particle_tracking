function  frequency = center_rotation(folder,En)
%% Use images and angle of rotation to find how the center of rotation moves
filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);

%% Get number of rotation steps and numbers of the files with the images. 
i_rot = find(diff(Rotation_step(1,:))~=0);
nb_rotation_steps = Rotation_step(1,i_rot);
n_files = length(nb_rotation_steps);
%% Saving variables
xot = zeros(1,n_files);
yot = zeros(1,n_files);
Rt = zeros(1,n_files);
%% Get initial range of xo, yo and R from first file
file_for_edge = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,Fn_imafile(i_rot(1)));
load(file_for_edge);
%get positions of one image
px = pxs(pxs(:,1) > 0);
py = pys(pys(:,1) > 0);
xx = repmat(px,1,length(px));
yy = repmat(py,1,length(py));
dx = xx - xx';
dy = yy - yy';
R2 = dx.^2+dy.^2;
[mR2, I] = max(R2(:));
[a, b] = ind2sub(size(R2),I);
xo_ini = min(px([a b])) + abs(dx(a,b))/2;
yo_ini = min(py([a b])) + abs(dy(a,b))/2;
R_ini = sqrt(mR2)/2;


%% Main loop over positions after each rotation to get center each time
for ii = 1:n_files
    file_number = Fn_imafile(i_rot(ii));
    [xot(ii), yot(ii), Rt(ii)] = get_circle_center(folder,En,file_number,xo_ini,yo_ini,R_ini);
end

%% Fit xot yot data to sin cos.
fit_x = fit(nb_rotation_steps',xot','Fourier1');
fit_y = fit(nb_rotation_steps',yot','Fourier1');

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

frequency = fit_x.w;

file_for_center = sprintf('%sCenter_%i.mat',filedirectory,En);
save(file_for_center,'xot','yot','Rt','fit_x','fit_y','phi_x','phi_y','A_x','A_y','frequency');