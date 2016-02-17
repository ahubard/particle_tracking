function  ntheta = center_rotation(folder,En)
%% Use images and angle of rotation to find how the center of rotation moves

% folder = 1;
% En = 104;
D = 10;

filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);
images_file = sprintf('%sTo_get_center_%i.mat',filedirectory,En);

%% Get number of rotation steps and numbers of the files with the images. 
original_image = 1;     % image to complete
angle_aperture = 2.35;
rot_angle = 0.0031914; %from minimizing distance btw centers in compare_rotation function.
i_final = find(diff(Rotation_step(1,:))~=0);
i_initial = find(diff(Rotation_step(1,:))~=0)+1;
i_final = i_final(original_image:end);
i_initial = i_initial(original_image:end);
nb_rotation_steps_btw_avalanches = Rotation_step(1,i_initial)-Rotation_step(1,i_final);%fix this? 
Total_rotation_steps = cumsum(nb_rotation_steps_btw_avalanches);
Dangle_btw_avalanches = nb_rotation_steps_btw_avalanches*rot_angle;
Tot_angle = cumsum(Dangle_btw_avalanches);
ntheta = length(Tot_angle);

%% Get first image or initial state before any rotations. 
before_fn = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',...
        folder,folder,En,En,Fn_imafile(i_final(1)));
load(before_fn,'IMA');
rot_images = zeros(size(IMA,1),size(IMA,2),ntheta);
rot_images(:,:,1) = sum(IMA,3)/sum(IMA(:));
for ii = 1:ntheta
after_fn = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',...
folder,folder,En,En,In_imafile(i_initial(ii)));
load(after_fn,'IMA')
rot_images(:,:,ii+1) = sum(IMA,3)/sum(IMA(:));
end
save(images_file,'Total_rotation_steps','Tot_angle','rot_images');