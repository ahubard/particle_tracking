%% Find the positions of the particles that dont appear in the picture. 

% Data keeping files
folder = 1; 
En = 103;


filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);
rotation_center_file = sprintf('%srotation_center_%i.mat',filedirectory,En);
load(rotation_center_file);

%% E 

angle_aperture = 2.35;
rot_angle = 0.0032;
i_final = find(diff(Rotation_step(1,:))~=0);
i_initial = find(diff(Rotation_step(1,:))~=0)+1;d_R = Rotation_step(1,i_final)-Rotation_step(1,i_initial);
th = d_R*rot_angle;


for ii = 1:length(i_initial)
    before_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,Fn_imafile(i_final(ii)));
    after_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,In_imafile(i_initial(ii)));
    load(before_fn,'pxs','pys','Npf');
    pxb(1:Npf(351),ii) = pxs(1:Npf(351),351);
    pyb(1:Npf(351),ii) = pys(1:Npf(351),351);
    load(after_fn,'pxs','pys','Npf')
    pxa(1:Npf(1),ii) = pxs(1:Npf(1),1);
    pya(1:Npf(1),ii) = pys(1:Npf(1),1);
end

%% Check to see if it matches.


for ii 