function [x_rot, y_rot] = rotation_composition(x,y,o_t,r_t,A_x,A_y,phi_x,phi_y, frequency,a0_x,a0_y)
% rotation_composition rotates coordinates x,y from rotation_image to the
% reference of the original image. Adds a correction over the fact that the
% rotation center of the image moves as the drum rotates.
%o_t is the rotaton step of the image to compose, r_t the one of the image
%to rotate, A_x, A_y, phi_x, phi_y, frequency,a0_x and a0_y are the
%function parameters of how the center move.


%% Get the number of rotation steps and recreate the center in each step.
n_rot = r_t - o_t;
t = [o_t r_t];  %List of rotation steps between o_t and r_t
%Positions of the center btw o_t and r_t
xo = A_x*sin(frequency*t+phi_x)+a0_x;   
yo = A_y*cos(frequency*t+phi_y)+a0_y;
%Get difference btw xo and yo to apply the rotation matrices
theta = n_rot*frequency; %Angle that will be rotated
[x_rot, y_rot] = rot_me(theta,x-xo(2),y-yo(2));
x_rot = x_rot + xo(1);
y_rot = y_rot + yo(1);



