filenumbers = [ 15 16 17 18 19 20 21 22 23 103 104 105 106 107 108 109 ]; %Files that contain the info
Nofiles = length(filenumbers);
nf = 10;
Nbparticles = 6250;
R = 578;
half_circle_cm = [0; -4*R/(3*pi)];
rcm2 = sum(half_circle_cm.^2);

rot_angle = .21*pi/180;
initialangle = 29*pi/180;

filename = sprintf('Avalanches_%i.mat',filenumbers(nf));



load(filename,'git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','DLength','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary','diff_Center_mass',...
    'Displacement_File_nb', 'Participation', 'In_imafile','Fn_imafile','in_trackedfile');

Nb_rot = diff(Rotation_step(1,:));
CM = zeros(2,2*(Number_Avalanches-1)-1);
rotationangle(1) = Initial_Angle(1)*pi/180+initialangle;
rotation_matrix = [cos(rotationangle(1)) -sin(rotationangle(1)) ; ...
sin(rotationangle(1)) cos(rotationangle(1))];
CM(:,1) = rotation_matrix*half_circle_cm;
jj = 2;
for ii = 2:Number_Avalanches-1
    CM(2,jj) = CM(2,jj-1)+(Dheight(ii)/Nbparticles);
    CM(1,jj) = CM(1,jj-1)+(DLength(ii)/Nbparticles);
    
    rotationangle(ii) = rot_angle*Nb_rot(ii);
    rotation_matrix = [cos(rotationangle(ii)) -sin(rotationangle(ii)) ; ...
    sin(rotationangle(ii)) cos(rotationangle(ii))];
    CM(:,jj+1) = rotation_matrix*CM(:,jj);
  
    jj = jj+2;
   
end

    
    






