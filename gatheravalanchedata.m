%% Gather data from series of experiments to do the statistics.

filenumbers = [14 104 15  106 17 18 109 19 20 23 103  105  107 16  ]; %Files that contain the info
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





%% Load and read files
for nf = 1:Nofiles
    filename = sprintf('Avalanches_%i.mat',filenumbers(nf));
    
    clear('git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step');
    %'spectrum_particles','spectrum_displacement','spectrum_energy','spectrum_potential',...
    
    
    load(filename,'git_version','Number_Avalanches','Noavalanches','Avalanche_time', ...
        'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
        'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
        'mat_particles','mat_displacement','mat_energy','mat_potential',...
        'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
        'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
        'Initial_Angle','Final_Angle','Rotation_step');
    %'spectrum_particles','spectrum_displacement','spectrum_energy','spectrum_potential',...
    
    Totalavalanches = Totalavalanches + Number_Avalanches;
    filenumber = [filenumber nf*ones(1,Number_Avalanches)];
    dSteps = [dSteps diff(Rotation_step(1,:)) -1];
    T = [T Avalanche_duration];
    
    Size_Particles = [Size_Particles Avalanche_particles];
    Size_Length_path = [Size_Length_path Avalanche_displacement];
    Size_Energy = [Size_Energy Avalanche_energy];
    Size_Potential = [Size_Potential Avalanche_potential];
    
    
    
    shape_Particles = [shape_Particles Normalized_particles];
    shape_Length_path = [shape_Length_path Normalized_avalanche];
    shape_Energy = [shape_Energy Normalized_energy];
    shape_Potential = [shape_Potential Normalized_potential];
    
    Particles_t = [Particles_t mat_particles];
    Length_path_t = [Length_path_t mat_displacement];
    Energy_t = [Energy_t mat_energy];
    Potential_t = [Potential_t mat_potential];
    
    Particles_Correlation = [Particles_Correlation correlation_particles];
    Length_path_Correlation = [Length_path_Correlation correlation_displacement];
    Energy_Correlation = [Energy_Correlation correlation_energy];
    Potential_Correlation = [Potential_Correlation correlation_potential];
    
    %     Power_spectrum_Particles = [Power_spectrum_Particles spectrum_particles];
    %     Power_spectrum_Length_path = [Power_spectrum_Length_path spectrum_displacement];
    %     Power_spectrum_Energy = [Power_spectrum_Energy spectrum_energy];
    %     Power_spectrum_Potential = [Power_spectrum_Potential spectrum_potential];
    
    Total_Particles = [Total_Particles NoParticles_moved];
    Total_Position_change = [Total_Position_change DELTAR];
    Itheta = [Itheta Initial_Angle];
    Ftheta = [Ftheta Final_Angle];
    Maximal_particle_displacement = [Maximal_particle_displacement Max_particle_dis];
    Total_Height_change = [Total_Height_change Dheight];
    
end
Delta_theta = Itheta-Ftheta;
ii = find(dSteps>-1);
T = T/fps;
Size_Length_path = Size_Length_path/D;
Maxd = Maxd/D;
TSd = TSd/D;
H = H/D;
U = H*g*d;
Energy = Energy/(D^2);
Energyu = Energy*((d/D)^2)*(fps^2) ;
Itheta = Itheta+alpha;
Ftheta = Ftheta+alpha;
ir = ones(size(dSteps));
ir(1) = 0;
inr = find(dSteps<0);
ir(inr) = 0;
ir(inr+1)= 0;
ir = find(ir);
