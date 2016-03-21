filenumbers = [ 15 16 17 18 19 20 21 22 23 103 104 105 106 107 108 109 ]; %Files that contain the info
Nofiles = length(filenumbers);
nf = 10;
initial = 1;
kind = 2; 
%% Load files
   En = filenumbers(nf);
    filedirectory = sprintf('/Users/Aline/Documents/Research/MATLAB/particle_tracking/');
    
    filename = sprintf('%sAvalanches_%i_%i.mat',filedirectory,En,kind);
    file_SCM = sprintf('%sSurface_CM_%i_%i.mat',filedirectory,En,kind);
    file_Potential = sprintf('%sPotential_Energy_%i_%i.mat',filedirectory,En,kind);
    file_CM = sprintf('%sCenter_of_Mass_%i.mat',filedirectory, En);
%filename = sprintf('Avalanches_%i.mat',filenumbers(nf));
load(filename,'git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','DLength','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary','diff_Center_mass',...
    'Displacement_File_nb', 'Participation', 'In_imafile','Fn_imafile','in_trackedfile');
%load(file_Potential);
load(file_CM);
if(En <100)

    R = 579.9;
else
 
    R = 617.7;
    
end
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile);

%% Get parameter of the rotation
[Amp, a0_x, a0_y, phi_x] = modify_center_func(1,En);
phi_y = phi_x +0.1;
%% Initial conditions
Nb = N_particles(1);
CM = zeros(2,2*(Number_Avalanches-initial)-1);
CM_from_angles = zeros(2,2*(Number_Avalanches-initial)-1);

rot_angle = .0032;
initialangle = 29*pi/180;
ini_CM = [x_cm(1) y_cm(1)];
t1 = Avalanche_time{1,initial};
[CM(1,1), CM(2,1)] = rot_me(-initialangle, ini_CM(1), ini_CM(2)); %FROM IMAGE
 %This is the center of mass of the half circle with the origin at a0
half_circle_cm = [0; 4*R/(3*pi)] +[a0_x; a0_y];

rcm2 = sum(half_circle_cm.^2);



ni = initial:Number_Avalanches;
th_ini = Initial_Angle(ni)*pi/180;
th_fin = Final_Angle(ni)*pi/180;
[CM_from_angles(1,1),CM_from_angles(2,1)] = ...
    rot_me(-th_ini(1), half_circle_cm(1),half_circle_cm(2));

Dheight = Dheight(ni);
DLength = DLength(ni);
Total_displacement = sqrt(Dheight.^2+DLength.^2);
rot_i = Rotation_step(1,ni);
Nb_rot = diff(Rotation_step(1,(ni)));
rotationangle = zeros(size(Nb_rot));
rotationangle(1) = 0;
%%
jj = 2;
for ii = 1:Number_Avalanches-initial-1
    
     CM(2,jj) = CM(2,jj-1)-(Dheight(ii)/Nb);
     CM(1,jj) = CM(1,jj-1)+(DLength(ii)/Nb); 
%      rotationangle(ii+1) = rot_angle*Nb_rot(ii);
%      rotation_matrix = [cos(rotationangle(ii+1)) -sin(rotationangle(ii+1)); ...
%      sin(rotationangle(ii+1)) cos(rotationangle(ii+1))];
%      CM(:,jj+1) = rotation_matrix*CM(:,jj);
    [CM(1,jj+1),CM(2,jj+1)] = ...
        rotation_composition(CM(1,jj),CM(2,jj),rot_i(ii),rot_i(ii+1),...
        Amp,Amp,phi_x,phi_y,-rot_angle,a0_x,a0_y);

     [CM_from_angles(1,jj),CM_from_angles(2,jj)] =...
         rot_me(-th_fin(ii),half_circle_cm(1),half_circle_cm(2));
     [CM_from_angles(1,jj+1),CM_from_angles(2,jj+1)] =...
         rot_me(-th_ini(ii+1),half_circle_cm(1),half_circle_cm(2));
     jj = jj+2;
   
end

x1 = cumsum(Total_displacement(1:end-1));
y1 = cumsum(rot_angle*Nb_rot*4*R/(3*pi));   
    






