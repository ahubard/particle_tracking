%% Gather data from series of experiments to do the statistics.

filenumbers = [ 15 16 17 18 19 20 21 22 23 103 104 105 106 107 108 109 ]; %Files that contain the info
Nofiles = length(filenumbers);
%% Parameters from the experiment
fps = 694.4444;
D = 10;
alpha = 29;
g = 9.8;
d = 12e-3;

%% Create variables
Totalavalanches = 0;
filenumber = [];
dSteps = [];
T = [];

Particles_t = [];
Length_path_t = [];
Energy_t = [];
Potential_t = [];

Particles_Correlation = [];
Length_path_Correlation = [];
Energy_Correlation = [];
Potential_Correlation = [];

Size_Particles = [];
Size_Length_path = [];
Size_Energy = [];
Size_Potential = [];

shape_Particles = []; 
shape_Length_path =[];
shape_Energy = [];
shape_Potential = [];

% Power_spectrum_Particles = [];
% Power_spectrum_Length_path = [];
% Power_spectrum_Energy = [];
% Power_spectrum_Potential = [];

Total_Particles = [];
Total_Position_change = [];
Itheta = [];
Ftheta = [];
Maximal_particle_displacement = [];
Total_Height_change = [];
Particles_Over_the_boundary = [];





%% Load and read files
for nf = 1:Nofiles
    filename = sprintf('Avalanches_%i.mat',filenumbers(nf));
    
    clear('git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary');
    %'spectrum_particles','spectrum_displacement','spectrum_energy','spectrum_potential',...
    
    
    load(filename,'git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary');
    %'spectrum_particles','spectrum_displacement','spectrum_energy','spectrum_potential',...
    
    ii = 2:Number_Avalanches;

    
    Particles_Over_the_boundary = [Particles_Over_the_boundary Nb_boundary(:,ii)];
    Totalavalanches = Totalavalanches + Number_Avalanches-1;
    filenumber = [filenumber nf*ones(1,Number_Avalanches-1)];
    dSteps = [dSteps diff(Rotation_step(1,:)) -1];
    T = [T Avalanche_duration(ii)];
    
    Size_Particles = [Size_Particles Avalanche_particles(ii)];
    Size_Length_path = [Size_Length_path Avalanche_displacement(ii)];
    Size_Energy = [Size_Energy Avalanche_energy(ii)];
    Size_Potential = [Size_Potential Avalanche_potential(ii)];
    
    
    
    shape_Particles = [shape_Particles Normalized_particles(:,ii)];
    shape_Length_path = [shape_Length_path Normalized_avalanche(:,ii)];
    shape_Energy = [shape_Energy Normalized_energy(:,ii)];
    shape_Potential = [shape_Potential Normalized_potential(:,ii)];
    
    Particles_t = [Particles_t mat_particles(:,ii)];
    Length_path_t = [Length_path_t mat_displacement(:,ii)];
    Energy_t = [Energy_t mat_energy(:,ii)];
    Potential_t = [Potential_t mat_potential(:,ii)];
    
    Particles_Correlation = [Particles_Correlation correlation_particles(:,ii)];
    Length_path_Correlation = [Length_path_Correlation correlation_displacement(:,ii)];
    Energy_Correlation = [Energy_Correlation correlation_energy(:,ii)];
    Potential_Correlation = [Potential_Correlation correlation_potential(:,ii)];
    
    %     Power_spectrum_Particles = [Power_spectrum_Particles spectrum_particles];
    %     Power_spectrum_Length_path = [Power_spectrum_Length_path spectrum_displacement];
    %     Power_spectrum_Energy = [Power_spectrum_Energy spectrum_energy];
    %     Power_spectrum_Potential = [Power_spectrum_Potential spectrum_potential];
    
    Total_Particles = [Total_Particles NoParticles_moved(ii)];
    Total_Position_change = [Total_Position_change DELTAR(ii)];
    Itheta = [Itheta Initial_Angle(ii)];
    Ftheta = [Ftheta Final_Angle(ii)];
    Maximal_particle_displacement = [Maximal_particle_displacement Max_particle_dis(ii)];
    Total_Height_change = [Total_Height_change Dheight(ii)];
    
end

 ii_bound_left = find(Particles_Over_the_boundary(1,:) > 0);
 ii_bound_right = find(Particles_Over_the_boundary(2,:) > 1);
 ii_boundary = find(((Particles_Over_the_boundary(2,:) > 1)+(Particles_Over_the_boundary(1,:) > 0)));
 ii_boundary_both = find(((Particles_Over_the_boundary(2,:) > 1).*(Particles_Over_the_boundary(1,:) > 0)));
 ii_internal_left = find(Particles_Over_the_boundary(1,:) <= 0);
 ii_internal_right = find(Particles_Over_the_boundary(2,:) <= 1);
 ii_none_boundary = find(((Particles_Over_the_boundary(1,:) <= 0).*(Particles_Over_the_boundary(2,:) <= 1)));
Delta_theta = Itheta-Ftheta;
ii = find(dSteps>-1);
T = T/fps;
% Size_Length_path = Size_Length_path/D;
% Maxd = Maxd/D;
% TSd = TSd/D;
% H = H/D;
% U = H*g*d;
% Energy = Energy/(D^2);
% Energyu = Energy*((d/D)^2)*(fps^2) ;
Itheta = Itheta+alpha;
Ftheta = Ftheta+alpha;
ir = ones(size(dSteps));
ir(1) = 0;
inr = find(dSteps<0);
ir(inr) = 0;
ir(inr+1)= 0;
ir = find(ir);
