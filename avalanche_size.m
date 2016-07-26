function Number_Avalanches = avalanche_size(folder,En,kind)
%% Measure size and shape of avalanches

%% Separate the avalanches according to kind criterium. 
% If kind = 0 the separation is the same as in displacement. Meaning all
% the movement that correspond to the same continuos movie are part of the
% same avalanche.
% If kind = 1, an avalanche comprises everything that happens after one
% rotation. In general less avalanches than in kind = 0;
% If kind = 2, an avalanche starts when particles move and ends when they
% stop for more than minstepsbtwav. In general more avalanches than in kind=0

%% Version Control
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');

%% Gather the basic info
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
%avanofile = sprintf('Avanonestep%i.mat',En);
image_file_length = 350;
if (exist(avanofile,'file'))
    load (avanofile,'count','avan','alpha','maxT','navfile');
    if (kind < 2)
        changefileindex = find(diff(avan(1,navfile))>1);
        initialfileindex = navfile([1 changefileindex(1:end-1)+1]);%first file of set
        finalfileindex = navfile(changefileindex);%last file of set.
        maxT = max(1+finalfileindex - initialfileindex)*350;
    end
    
else
    save(sprintf('AWarning. The file: %s does not exist.mat',avanofile));
    error('Error, avanofile does not exist');
end

alpha = alpha*pi/180;

if(En < 100)

    R = 579.9;
else
 
    R = 617.7;
    
end

[Amp, xo, yo, phi_x] = modify_center_func(folder,En);
phi_y = phi_x +0.1;



D = 10;
Ly = 400;
yo = Ly-yo; 
%% Initialize variables.
mat_displacement = zeros(maxT,4*count);
mat_particles = zeros(maxT,4*count);
mat_energy = zeros(maxT,4*count);
mat_potential = zeros(maxT,4*count);

correlation_particles = zeros(2*maxT-1,4*count);
correlation_displacement = zeros(2*maxT-1,4*count);
correlation_energy  = zeros(2*maxT-1,4*count);
correlation_potential = zeros(2*maxT-1,4*count);

Avalanche_displacement = zeros(1,4*count);
Avalanche_energy = zeros(1,4*count);
Avalanche_duration = zeros(1,4*count);
Avalanche_particles = zeros(1,4*count);
Avalanche_potential = zeros(1,4*count);

Normalized_avalanche = zeros(101,4*count);
Normalized_energy = zeros(101,4*count);
Normalized_particles = zeros(101,4*count);
Normalized_potential = zeros(101,4*count);

Nb_bins = ceil(2*R/(D));
diff_Center_mass = zeros(Nb_bins,count);

diff_CMass_t = cell(1,count);
Total_Potential_Energy = cell(1,count);
% spectrum_displacement = zeros(maxT/2,4*count);
% spectrum_particles = zeros(maxT/2,4*count);
% spectrum_energy = zeros(maxT/2,4*count);
% spectrum_potential = zeros(maxT/2,4*count);
Nb_boundary = zeros(2,count);
Number_Avalanches = 0;
Noavalanches = zeros(1,count);

%Total x and y change in camera coordinate system
Delta_x = zeros(1,4*count);
Delta_y = zeros(1,4*count);

DELTAR = zeros(1,4*count);
Dheight = zeros(1,4*count);
DLength = zeros(1,4*count);
NoParticles_moved = zeros(1,4*count);
Max_particle_dis = zeros(1,4*count);
%Initial_CM = zeros(2,4*count);
%Final_CM = zeros(2,4*count);
%Num_part_ini = zeros(1,4*count);
%Num_part_end = zeros(1,4*count);
Initial_Angle = zeros(1,4*count);
Final_Angle = zeros(1,4*count);
Rotation_step = zeros(2,4*count);
Avalanche_time = cell(2,count);

%% Pass of the information from displacement files for later analisys
Nb_Tracked = zeros(1,count);
Participation = zeros(1,count);
In_imafile = zeros(1,count);
Fn_imafile = zeros(1,count);
in_trackedfile = zeros(1,count);
Displacement_File_nb = zeros(1,count);
%% Cutoff values

cutoffperparticle = 0.01;    %If drfil (filtered particle displacement) is smaller
% cutof then the particle didn move.
cutoffofsum = 0.06;
windowSize = 5;
minstepsbtwav = 10; %same as increment. Minimal amount of steps whit no displacement to say the avalanche is over.
b = (1/minstepsbtwav)*ones(1,minstepsbtwav);
a = 1;


%% Main loop

for nf = 1:count-1
    fnn =sprintf('%sDisplacement_%i.mat',filedirectory,nf);
    
    clear('windowSize', 'PX','PY','Nb_over_boundary','dh','dl','drraw2','drfil2','diskmove','increment','participationratio','NEn','initial','final');
  
    if(exist(fnn,'file'))
        load(fnn,'windowSize', 'PX','PY','Nb_over_boundary','dh','dl','drraw2','drfil2','diskmove','increment','participationratio','NEn','initial','final');
    else
        save(sprintf('Warning. The file: %s does not exist.mat',fnn));
        error('Error, Displacementfile does not exist');
    end
    
    if (diskmove)%Check there was no error in the file.
        
        drmask = (drfil2>=cutoffperparticle);  %Use filter dr to determine which particles moved at each time step.
        particles = sum(drmask);
        energy_avalanche = sum(drraw2(:,ceil(windowSize/2):size(drraw2,2)-floor(windowSize/2)).*drmask);
        vel_avalanche = sum(sqrt(drraw2(:,ceil(windowSize/2):size(drraw2,2)-floor(windowSize/2))).*drmask);
        dif_potential = sum(dh(:,ceil(windowSize/2):size(dh,2)-floor(windowSize/2)).*drmask);
        findavalanche = filter(b,a,(energy_avalanche>cutoffofsum)); %use both numbers to decide, if they are moving they should move for at least 10 time steps, if they stop the same
        % Find how many avalancheso13 per file findavalanche = 0 means no
        % avalanche there.
        
        %Posible beginings and endings of avalanches
        t1 = find(diff(findavalanche>0)==1);  %Time steps avalanche begins
        if(isempty(t1))
            t1 = 1;
        end
        
        t2 = min(find(diff(findavalanche>0)==-1)+1,length(energy_avalanche)); %Time steps avalanche ends.
        if(isempty(t2))
            t2 = length(energy_avalanche);
        end
        if (t1(1)>= t2(1)) %In case avalanche starts at first step.
            t1 = [1 t1];
        end
        
        if(t2(end)<=t1(end))%In case avalanche finish at the last frame of the file.
            t2 = [t2 length(energy_avalanche)];
        end
        
        % Only one avalanche per file in kind =0,1
        if( kind < 2)
            t1 = t1(1);
            t2 = t2(end);
        end
        
        nb_avalanches = length(t1);
        ii_time = zeros(1,nb_avalanches);
        
        [~, ii_particles_for_potential] = sort(PY(:,1));
        
        for na = 1:nb_avalanches
            if(sum((findavalanche(t1(na):min(t2(na),length(energy_avalanche))))>=1-eps)) %Check if there is indeed avalanches between t1-t2
                
                ii_time(na) = 1;

                
                Number_Avalanches = Number_Avalanches+1;
                %General data of file
                Nb_Tracked(Number_Avalanches) = size(PX,1);
                Displacement_File_nb(Number_Avalanches) = nf;
                Participation(Number_Avalanches) = participationratio;
                In_imafile(Number_Avalanches) = initial;
                Fn_imafile(Number_Avalanches) = final;
                in_trackedfile(Number_Avalanches) = NEn;                
                
               
                Rotation_step(:,Number_Avalanches) = avan(2,[initial; final]);
               

                
                %get the relevant data
                particlesmoving_t = particles(t1(na):t2(na));
                displacement_t = vel_avalanche(t1(na):t2(na));
                energy_t = energy_avalanche(t1(na):t2(na));
                potential_t = dif_potential(t1(na):t2(na));
                
                %center of mass
%                 file_ava_begins = floor(t1(na)/image_file_length)+initial;
%                 time_in_file_1 = mod(t1(na) - 1,image_file_length) + 1;
%                 [x_whole, y_whole] = create_bottom(folder,En,file_ava_begins,time_in_file_1);
%                 Initial_CM(:,Number_Avalanches) = [mean(x_whole) mean(y_whole)];
%                 Num_part_ini(Number_Avalanches) = length(x_whole);
%                 
%                 file_ava_end = floor(t2(na)/image_file_length)+initial;
%                 time_in_file_2 = mod(t2(na) - 1,image_file_length) + 1;
%                 [x_whole, y_whole] = create_bottom(folder,En,file_ava_end,time_in_file_2);
%                 Final_CM(:,Number_Avalanches) = [mean(x_whole); mean(y_whole)];
%                 Num_part_end(Number_Avalanches) = length(x_whole);
                
                
                %Sizes and duration
                Avalanche_particles(Number_Avalanches) = sum(particlesmoving_t);
                Avalanche_displacement(Number_Avalanches) = sum(displacement_t);
                Avalanche_energy(Number_Avalanches) = sum (energy_t);
                Avalanche_potential(Number_Avalanches) = sum(potential_t);
                
                deltat = t2(na)-t1(na);
                Avalanche_duration(Number_Avalanches) = deltat+2; %Adding the frame before and the frame after for completion.
                
                %Find normalized avalanches
                particlesnormalized = interp1(([ 0 windowSize+(0:deltat) deltat+2*windowSize])/(deltat+2*windowSize),[0 particlesmoving_t 0],(0:.01:1),'pchip');
                Normalized_particles(:,Number_Avalanches) = particlesnormalized;
                
                avalanchenormalized = interp1(([ 0 windowSize+(0:deltat) deltat+2*windowSize])/(deltat+2*windowSize),[0 displacement_t 0],(0:.01:1),'pchip');
                Normalized_avalanche(:,Number_Avalanches) = avalanchenormalized;
                
                energynormalized = interp1(([ 0 windowSize+(0:deltat) deltat+2*windowSize])/(deltat+2*windowSize),[0 energy_t 0],(0:.01:1),'pchip');
                Normalized_energy(:,Number_Avalanches) = energynormalized;
                
                potentialnormalized = interp1(([ 0 windowSize+(0:deltat) deltat+2*windowSize])/(deltat+2*windowSize),[0 potential_t 0],(0:.01:1),'pchip');
                Normalized_potential(:,Number_Avalanches) = potentialnormalized;
                
                
                %Fill data with zeros to get power spectrum
                particlesmoving_t(deltat-1:maxT) = 0;
                displacement_t(deltat-1:maxT) = 0;
                energy_t(deltat-1:maxT) = 0;
                potential_t(deltat-1:maxT) = 0;
                
                mat_particles(:,Number_Avalanches) = particlesmoving_t;
                mat_displacement(:,Number_Avalanches) = displacement_t;
                mat_energy(:,Number_Avalanches) = energy_t;
                mat_potential(:,Number_Avalanches) = potential_t;
                
                %correlations
                correlation_particles(:,Number_Avalanches) = xcorr(particlesmoving_t);
                correlation_displacement(:,Number_Avalanches) = xcorr(displacement_t);
                correlation_energy(:,Number_Avalanches) = xcorr(energy_t);
                correlation_potential(:,Number_Avalanches) = xcorr(potential_t);
                
                %
                %                 % Fourier transform data
                %                 Fparticlesmoving_t = abs(fft(particlesmoving_t))/(maxT/2);
                %                 Fdisplacement_t = abs(fft(displacement_t))/(maxT/2);
                %                 Fenergy_t = abs(fft(energy_t))/(maxT/2);
                %                 Fpotential_t = abs(fft(potential_t))/(maxT/2);
                %                 %PowerSpectrum
                %                 spectrum_particles(:,Number_Avalanches) = Fparticlesmoving_t(1:maxT/2).^2;
                %                 spectrum_displacement(:,Number_Avalanches) = Fdisplacement_t(1:maxT/2).^2;
                %                 spectrum_energy(:,Number_Avalanches) = Fenergy_t(1:maxT/2).^2;
                %                 spectrum_potential(:,Number_Avalanches) = Fpotential_t(1:maxT/2).^2;
                
                %Final-Initial data
                particlesthatmoved = (((PX(diskmove,t2(na))-PX(diskmove,t1(na))).^2+...
                    (PY(diskmove,t2(na))-PY(diskmove,t1(na))).^2) >1); % More than one pixel
                
                DELTAR(Number_Avalanches) = sum(sqrt(((PX(diskmove,t2(na))-PX(diskmove,t1(na))).^2+...
                    (PY(diskmove,t2(na))-PY(diskmove,t1(na))).^2).*particlesthatmoved));
                
                Delta_x(Number_Avalanches) = sum((PX(diskmove,t2(na)) - PX(diskmove,t1(na))).*particlesthatmoved);
                Delta_y(Number_Avalanches) = sum((PY(diskmove,t2(na)) - PY(diskmove,t1(na))).*particlesthatmoved);
                
                [x_before_aval, y_before_aval] = rot_me(alpha,PX(diskmove,t1(na)),PY(diskmove,t1(na)));
                [x_after_aval, y_after_aval] = rot_me(alpha,PX(diskmove,t2(na)),PY(diskmove,t2(na)));
                Dheight(Number_Avalanches) = sum((y_after_aval - y_before_aval).*particlesthatmoved);
                DLength(Number_Avalanches) = sum((x_after_aval - x_before_aval).*particlesthatmoved);
                
                NoParticles_moved(Number_Avalanches) = sum(particlesthatmoved);
                Nb_boundary(:,Number_Avalanches) = Nb_over_boundary;
                
                Max_particle_dis(Number_Avalanches) = max(sqrt((PX(diskmove,t2(na))-PX(diskmove,t1(na))).^2+...
                    (PY(diskmove,t2(na))-PY(diskmove,t1(na))).^2));
                
                [thetai, ~, xy_massi, ii_surfacei, bindexi] = estimate_angle(PX(:,t1(na)),PY(:,t1(na)),D,yo,Nb_bins, 1);
                Initial_Angle(Number_Avalanches) =  thetai;
                [thetaf, ~, xy_massf, ii_surfacef, bindexf] = estimate_angle(PX(:,t2(na)),PY(:,t2(na)),D,yo,Nb_bins,1);
                Final_Angle(Number_Avalanches) = thetaf;
                dummy_aux1 = zeros(Nb_bins, deltat+1);
                dummy_aux2 = zeros(1,deltat+1);
                
                for t = t1(na):t2(na)
                    dummy_aux1(:,t-t1(na)+1) = compare_C_Mass(PX(ii_surfacei,t),PY(ii_surfacei,t),bindexi,xy_massi,Nb_bins);
                    [~, ay] = rot_me(alpha,PX(ii_particles_for_potential ,t), PY(ii_particles_for_potential ,t));
                    dummy_aux2(t-t1(na)+1) = sum(ay);
                end
                
                diff_CMass_t{Number_Avalanches} = dummy_aux1;
                Total_Potential_Energy{Number_Avalanches} = dummy_aux2;
                aux_i = dummy_aux1(:,deltat+1);
                aux_f = compare_C_Mass(PX(ii_surfacef,t2(na)),PY(ii_surfacef,t2(na)),bindexf,xy_massf,Nb_bins);
                diff_Center_mass(:,Number_Avalanches) = aux_i+aux_f';
            end
            
            
        end
        Avalanche_time{1,nf} = t1(ii_time > 0);
        Avalanche_time{2,nf} = t2(ii_time > 0);
        Noavalanches(nf) = Number_Avalanches;
    end
    %fprintf('In file %i the number of avalanches is %i\n',nf,Number_Avalanches);
    
end

    
%% Resize matrices keeping only the data the has a "big" participation ratio
ikeep = find(Participation > 0.0025);

Number_Avalanches = length(ikeep);
Avalanche_duration = Avalanche_duration(ikeep);

Avalanche_particles = Avalanche_particles(ikeep);
Avalanche_displacement = Avalanche_displacement(ikeep);
Avalanche_energy = Avalanche_energy(ikeep);
Avalanche_potential = Avalanche_potential(ikeep);

Normalized_particles = Normalized_particles(:,ikeep);
Normalized_avalanche = Normalized_avalanche(:,ikeep);
Normalized_energy = Normalized_energy(:,ikeep);
Normalized_potential = Normalized_potential(:,ikeep);

mat_particles = mat_particles(:,ikeep);
mat_displacement = mat_displacement(:,ikeep);
mat_energy = mat_energy(:,ikeep);
mat_potential = mat_potential(:,ikeep);

%correlations
correlation_particles = correlation_particles(:,ikeep);
correlation_displacement = correlation_displacement(:,ikeep);
correlation_energy = correlation_energy(:,ikeep);
correlation_potential = correlation_potential(:,ikeep);
diff_Center_mass = diff_Center_mass(:,ikeep);
diff_CMass_t = diff_CMass_t(1,ikeep);

Total_Potential_Energy = Total_Potential_Energy(ikeep);

% spectrum_particles = spectrum_particles(:,ikeep);
% spectrum_displacement = spectrum_displacement(:,ikeep);
% spectrum_energy = spectrum_energy(:,ikeep);
% spectrum_potential = spectrum_potential(:,ikeep);
Nb_boundary = Nb_boundary(:,ikeep);
DELTAR = DELTAR(ikeep);
Delta_x = Delta_x(ikeep);
Delta_y = Delta_y(ikeep);

Dheight = Dheight(ikeep);
DLength = DLength(ikeep);
NoParticles_moved = NoParticles_moved(ikeep);
Max_particle_dis = Max_particle_dis(ikeep);

% Initial_CM = Initial_CM(:,ikeep);
% Final_CM =Final_CM(:,ikeep);
% Num_part_ini = Num_part_ini(ikeep);
% Num_part_end = Num_part_end(ikeep);
Initial_Angle = Initial_Angle(ikeep);
Final_Angle = Final_Angle(ikeep);
Rotation_step = Rotation_step(:,ikeep);

Nb_Tracked = Nb_Tracked(ikeep);
Displacement_File_nb = Displacement_File_nb(ikeep);
Participation = Participation(ikeep);
In_imafile = In_imafile(ikeep);
Fn_imafile = Fn_imafile(ikeep);
in_trackedfile = in_trackedfile(ikeep);
%% Save results to file
file_save = sprintf('%sAvalanche_size_%i_%i.mat',filedirectory,En,kind);
file_CM = sprintf('%sSurface_CM_%i_%i.mat',filedirectory,En,kind);
file_Potential = sprintf('%sPotential_Energy_%i_%i.mat',filedirectory,En,kind);

% save(file_save,'git_version','maxT','Number_Avalanches','Noavalanches','Avalanche_time', ...
%     'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
%     'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
%     'mat_particles','mat_displacement','mat_energy','mat_potential',...
%     'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
%     'DELTAR','Dheight','DLength','NoParticles_moved','Max_particle_dis',...
%     'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary','diff_Center_mass',...
%     'Displacement_File_nb', 'Participation', 'In_imafile','Fn_imafile','in_trackedfile',...
%     'Delta_x','Delta_y','Initial_CM','Final_CM','Num_part_ini','Num_part_end'); 
save(file_save,'git_version','maxT','Number_Avalanches','Noavalanches','Avalanche_time', ...
    'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
    'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
    'mat_particles','mat_displacement','mat_energy','mat_potential',...
    'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
    'DELTAR','Dheight','DLength','NoParticles_moved','Max_particle_dis',...
    'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary','diff_Center_mass',...
    'Displacement_File_nb', 'Participation', 'In_imafile','Fn_imafile','in_trackedfile',...
    'Delta_x','Delta_y'); 
% save(file_save,'git_version','MaxT','Number_Avalanches','Noavalanches','Avalanche_time', ...
%     'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
%     'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
%     'mat_particles','mat_displacement','mat_energy','mat_potential',...
%     'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
%     'spectrum_particles','spectrum_displacement','spectrum_energy','spectrum_potential',...
%     'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
%     'Initial_Angle','Final_Angle','Rotation_step');

save(file_CM,'git_version','diff_CMass_t','Avalanche_time','Displacement_File_nb','Delta_x','Delta_y');
%save(file_Potential,'git_version','Total_Potential_Energy','Initial_Angle',...
    'Final_Angle','Rotation_step','Dheight','DLength','Nb_Tracked',...
    'Delta_x','Delta_y','Initial_CM','Final_CM','Num_part_ini','Num_part_end');
