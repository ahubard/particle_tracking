filenumbers = [ 15 16 17 18 19 20 21 22 23 103 104 105 106 107 108 109 ]; %Files that contain the info
Nofiles = length(filenumbers);
nf = 2;
%12  7967
Nbparticles = 7500;
%Nbparticles = 8049;
%Nbparticles = 1;


filename = sprintf('Avalanches_%i.mat',filenumbers(nf));
avafile = sprintf('Avanonestep%02d.mat',filenumbers(nf));
load(avafile,'R');
%R = 578;
half_circle_cm = [0; -4*R/(3*pi)];
rcm2 = sum(half_circle_cm.^2);

rot_angle = .18*pi/180;
initialangle = 29*pi/180;


load(filename,'git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','DLength','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary','diff_Center_mass',...
    'Displacement_File_nb', 'Participation', 'In_imafile','Fn_imafile','in_trackedfile');
initial =1;
ni = initial:Number_Avalanches;
th_ini = Initial_Angle(ni)*pi/180;
Dheight = Dheight(ni);
DLength = DLength(ni);
th_fin = Final_Angle(ni)*pi/180;
Total_displacement = sqrt(Dheight.^2+DLength.^2);
Nb_rot = diff(Rotation_step(1,(ni)));
CM = zeros(2,2*(Number_Avalanches-initial)-1);
CM_from_angles = zeros(2,2*(Number_Avalanches-initial)-1);
rotationangle = zeros(size(Nb_rot));
rotationangle(1) = th_ini(1)+initialangle;
rotation_matrix = [cos(rotationangle(1)) -sin(rotationangle(1)) ; ...
sin(rotationangle(1)) cos(rotationangle(1))];
CM(:,1) = rotation_matrix*half_circle_cm;
CM_from_angles(:,1) = CM(:,1);
jj = 2;
for ii = 1:Number_Avalanches-initial-1
     
     CM(2,jj) = CM(2,jj-1)+(Dheight(ii)/Nbparticles);
     CM(1,jj) = CM(1,jj-1)+(DLength(ii)/Nbparticles); 
     rotationangle(ii+1) = rot_angle*Nb_rot(ii);
     rotation_matrix = [cos(rotationangle(ii+1)) -sin(rotationangle(ii+1)); ...
     sin(rotationangle(ii+1)) cos(rotationangle(ii+1))];
     CM(:,jj+1) = rotation_matrix*CM(:,jj);
     
     rot_matrix = [cos(th_fin(ii)+initialangle) -sin(th_fin(ii)+initialangle); ...
     sin(th_fin(ii)+initialangle) cos(th_fin(ii)+initialangle)];
     CM_from_angles(:,jj) = rot_matrix*half_circle_cm;
     rot_matrix = [cos(th_ini(ii+1)+initialangle) -sin(th_ini(ii+1)+initialangle); ...
     sin(th_ini(ii+1)+initialangle) cos(th_ini(ii+1)+initialangle)];
     CM_from_angles(:,jj+1) = rot_matrix*half_circle_cm;
 
    
     
         

     jj = jj+2;
   
end

 x1 = cumsum(Total_displacement(1:end-1));
y1 = cumsum(rot_angle*Nb_rot*4*R/(3*pi));   
    






