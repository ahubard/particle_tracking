%% Load data
filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
kind = 0;
alpha = 29*pi/180;

avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile);

file_CM = sprintf('%sCenter_of_Mass_%i.mat',filedirectory, En);
save(file_CM,'N_particles','x_cm','y_cm','N_rot_particles',...
    'x_from_rot_cm','y_from_rot_cm','git_version');


file_avalanches = sprintf('%sAvalanches_%i_%i.mat',filedirectory,En,kind);
load(file_avalanches);

file_center = sprintf('%sCenter_%i.mat',filedirectory,En);
load(file_center);

changefileindex = find(diff(avan(1,navfile))>1);  
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);%first file of set
finalfileindex = navfile(changefileindex);%last file of set. 

[~, ir] = intersect(initialfileindex,In_imafile);
x_cm = x_cm(ir);
y_cm = y_cm (ir);
N_rot_particles = N_rot_particles(ir);
%y_cm = 800-y_cm; %to be in the natural cartesian system of the camera before rotation

[x_cm, y_cm] = rot_me(-alpha,x_cm,y_cm);

pot_energy = y_cm * mean(N_particles);

dr = diff(Rotation_step(1,:));
dtheta = dr*alpha;




