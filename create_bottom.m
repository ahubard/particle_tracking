%% Find the positions of the particles that dont appear in the picture.
tic;
% Data keeping files
folder = 1;
En = 103;
D = 10;

filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);
rotation_center_file = sprintf('%srotation_center_%i.mat',filedirectory,En);
load(rotation_center_file);

%% E

angle_aperture = 2.35;
rot_angle = 0.0031914; %from minimizing distance btw centers in compare_rotation function.
i_final = find(diff(Rotation_step(1,:))~=0);
i_initial = find(diff(Rotation_step(1,:))~=0)+1;
nb_rotation_steps_btw_avalanches = Rotation_step(1,i_initial)-Rotation_step(1,i_final);%fix this? 
Dangle_btw_avalanches = nb_rotation_steps_btw_avalanches*rot_angle;
Tot_angle = cumsum(Dangle_btw_avalanches);

%% Use positions of files after rotation to find the particles before.
%Number of avalanches after rotating at least the aperture angle
n_rot = find(Tot_angle > angle_aperture, 1,'first'); 

%Load initial positions
before_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,Fn_imafile(i_final(1)));
load(before_fn,'pxs','pys','Npf');
pxo = pxs(1:Npf(351),351);%before rotation, after the avalanche finishes.
pyo = pys(1:Npf(351),351);
pxb = zeros(Npf(351),n_rot);
pyb = zeros(Npf(351),n_rot);
pxa = zeros(Npf(351),n_rot);
pya = zeros(Npf(351),n_rot);
    
for ii = 1:n_rot
    before_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,Fn_imafile(i_final(ii)));
    load(before_fn,'pxs','pys','Npf');
    pxb(1:Npf(351),ii) = pxs(1:Npf(351),351);%before rotation, after the avalanche finishes.
    pyb(1:Npf(351),ii) = pys(1:Npf(351),351);
    after_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,In_imafile(i_initial(ii))); 
    load(after_fn,'pxs','pys','Npf')
    pxa(1:Npf(1),ii) = pxs(1:Npf(1),1);%after rotation, before the avalanche starts.
    pya(1:Npf(1),ii) = pys(1:Npf(1),1);
end

%%
ii = n_rot;
%[~, trivialbondt1,trivialbondt2] = adjacent(pxa(:,ii-1),pya(:,ii-1),pxb(:,ii),pyb(:,ii),1);
%x = (pxa(trivialbondt1,ii-1) + pxb(trivialbondt2,ii))/2;
%y = (pya(trivialbondt1,ii-1) + pyb(trivialbondt2,ii))/2;
x = pxa(:,ii);
y = pya(:,ii);
yh = 400-R*sin(Dangle_btw_avalanches(ii));   %Minimal height of column to fill
i_zone = find(x > xo & y > (yh-D));
alpha = Tot_angle(ii+1);
rot_x = (x(i_zone)-xo)*cos(alpha)-(y(i_zone)-yo).*sin(alpha)+xo;
rot_y = (x(i_zone)-xo).*sin(alpha)+(y(i_zone)-yo).*cos(alpha)+yo;
izone = find(rot_y > (400-3*D));
x_old = rot_x(izone);
y_old = rot_y(izone);
plot(pxo,pyo,'.',x_old,y_old,'.');axis('equal');axis('ij');
drawnow;
for ii = n_rot-1:-1:2

%Get positions from the ii rotation
% [~, trivialbondt1,trivialbondt2] = adjacent(pxa(:,ii-1),pya(:,ii-1),pxb(:,ii),pyb(:,ii),1);
% x = (pxa(trivialbondt1,ii-1) + pxb(trivialbondt2,ii))/2;
% y = (pya(trivialbondt1,ii-1) + pyb(trivialbondt2,ii))/2;
x = pxa(:,ii);
y = pya(:,ii);
yh = 400-R*sin(Dangle_btw_avalanches(ii));   %Minimal height of column to fill

i_zone = find(x > xo & y > (yh-2*D)); % Keep only the one unafected by the avalanche.
alpha = Tot_angle(ii+1);
%rotate
rot_x = (x(i_zone)-xo)*cos(alpha)-(y(i_zone)-yo).*sin(alpha)+xo;
rot_y = (x(i_zone)-xo).*sin(alpha)+(y(i_zone)-yo).*cos(alpha)+yo;

izone = find(rot_y > (400-2*D)); %Keep only the one outside the original picture.
rot_x = rot_x(izone);
rot_y = rot_y(izone);

% Find the ones that overlap or come from both rotations
[~, trivialbondt1,trivialbondt2] = adjacent(x_old,y_old,tran_pxb,tran_pyb,2);
x_both = tran_pxb(trivialbondt2);
y_both = tran_pyb(trivialbondt2);
%Find the particles that are not in both sets. 
[x_old_e ix] = setdiff(x_old,x_old(trivialbondt1));
y_old_e = y_old(ix);
[tran_pxb_e ix] = setdiff(tran_pxb,tran_pxb(trivialbondt2));
tran_pyb_e = tran_pyb(ix);
%Positions from all previous rotations
x_old = [x_both; x_old_e; tran_pxb_e];
y_old = [y_both; y_old_e; tran_pyb_e];
end

toc