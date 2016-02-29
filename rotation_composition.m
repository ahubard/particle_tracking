function [x_rot, y_rot] = rotation_composition(x,y,original_image,rotation_image)
% rotation_composition rotates coordinates x,y from rotation_image to the
% reference of the original image. Adds a correction over the fact that the
% rotation center of the image moves as the drum rotates.

%% Load center Info
filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_for_center = sprintf('%sCenter_%i.mat',filedirectory,En);
load(file_for_center);
%Rotationstep of image to construct
o_t = nb_rotation_steps(original_image);
%Rotationstep of image that will be used to find the positions at o_t
r_t = nb_rotation_steps(rotation_image); 
n_rot = r_t-o_t;
t = o_t:sign(n_rot)*1:r_t;  %List of rotation steps between o_t and r_t
%Positions of the center btw o_t and r_t
xo = A_x*sin(frequency*t+phi_x)+fit_x.a0;   
yo = A_y*cos(frequency*t+phi_x-.2)+fit_y.a0;
%Get difference btw xo and yo to apply the rotation matrices
dxo = diff(xo);
dyo = diff(yo);
theta = (t(end-1:-1:1)-o_t)*frequency; %Angle that will be rotated
[rot_correction_x, rot_correction_y] = rot_me(theta,-dxo,-dyo);
rot_correction_x = sum(rot_correction_x) + xo(end);
rot_correction_y = sum(rot_correction_y) + yo(end);
[x_rot, y_rot]= rot_me(n_rot*frequency,x,y);
x_rot = x_rot + rot_correction_x;
y_rot = y_rot + rot_correction_y;



