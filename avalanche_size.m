function Number_Avalanches = avalanche_size(folder,En)
%% Gather the basic info

avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
%avanofile = sprintf('Avanonestep%i.mat',En);

if (exist(avanofile,'file'))
    load (avanofile,'count','avan','alpha','maxT');
else
    save(sprintf('AWarning. The file: %s does not exist.mat',avanofile));
    error('Error, avanofile does not exist');
end

[git_version, ~] = evalc('system(''git describe --dirty --alway'')');

salpha = sin(alpha*pi/180);
calpha = cos(alpha*pi/180);

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



% spectrum_displacement = zeros(maxT/2,4*count);
% spectrum_particles = zeros(maxT/2,4*count);
% spectrum_energy = zeros(maxT/2,4*count);
% spectrum_potential = zeros(maxT/2,4*count);
Nb_boundary = zeros(2,count);
Number_Avalanches = 0;
Noavalanches = zeros(1,count);
DELTAR = zeros(1,4*count);
Dheight = zeros(1,4*count);
NoParticles_moved = zeros(1,4*count);
Max_particle_dis = zeros(1,4*count);
Initial_Angle = zeros(1,4*count);
Final_Angle = zeros(1,4*count);
Rotation_step = zeros(2,4*count);
Avalanche_time = cell(2,count);

%% Pass of the information from displacement files for later analizys
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
    fnn =sprintf('/aline%i/rotdrum%i/o%02d/Displacement_%i.mat',folder,folder,En,nf);
    
    clear('windowSize', 'PX','PY','Nb_over_boundary','dh','drraw2','drfil2','diskmove','increment','participationratio','NEn','initial','final');
  
    if(exist(fnn,'file'))
        load(fnn,'windowSize', 'PX','PY','Nb_over_boundary','dh','drraw2','drfil2','diskmove','increment','participationratio','NEn','initial','final');
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
        findavalanche = filter(b,a,(energy_avalanche>cutoffofsum)); %use both numbers tdecide, if they are moving they should move for at least 10 time steps, if they stop the same
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
        
        
        for na = 1:length(t1)
            if( sum((findavalanche(t1(na):min(t2(na),length(energy_avalanche))))>=1-eps)) %Check if there is indeed avalanches between t1-t2
                
                %General data of file
                Displacement_File_nb(Number_Avalanches) = nf;
                Participation(Number_Avalanches) = participationratio;
                In_imafile(Number_Avalanches) = initial;
                Fn_imafile(Number_Avalanches) = final;
                in_trackedfile(Number_Avalanches) = NEn;                
                
                Number_Avalanches = Number_Avalanches+1;
                Rotation_step(:,Number_Avalanches) = avan(2,[initial; final]);
               

                
                %get the relevantdata
                particlesmoving_t = particles(t1(na):t2(na));
                displacement_t = vel_avalanche(t1(na):t2(na));
                energy_t = energy_avalanche(t1(na):t2(na));
                potential_t = dif_potential(t1(na):t2(na));
                
                
                
                %Sizes and duration
                Avalanche_particles(Number_Avalanches) = sum(particlesmoving_t);
                Avalanche_displacement(Number_Avalanches) = sum(displacement_t);
                Avalanche_energy(Number_Avalanches) = sum (energy_t);
                Avalanche_potential(Number_Avalanches) = sum(potential_t);
                
                deltat = t2(na)-t1(na);
                Avalanche_duration(Number_Avalanches) = deltat+3; %Adding the frame before and the frame after for completion.
                
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
                    (PY(diskmove,t2(na))-PY(diskmove,t1(na))).^2)>1); % More than one pixel
                
                DELTAR(Number_Avalanches) = sum(sqrt(((PX(diskmove,t2(na))-PX(diskmove,t1(na))).^2+...
                    (PY(diskmove,t2(na))-PY(diskmove,t1(na))).^2).*particlesthatmoved));
                
                Dheight(Number_Avalanches) = sum(((PY(diskmove,t2(na))*calpha+PX(diskmove,t2(na))*salpha)-...
                    (PY(diskmove,t1(na))*calpha+PX(diskmove,t1(na))*salpha)).*particlesthatmoved);
                
                
                NoParticles_moved(Number_Avalanches) = sum(particlesthatmoved);
                Nb_boundary(:,Number_Avalanches) = Nb_over_boundary;
                
                Max_particle_dis(Number_Avalanches) = max(sqrt((PX(diskmove,t2(na))-PX(diskmove,t1(na))).^2+...
                    (PY(diskmove,t2(na))-PY(diskmove,t1(na))).^2));
                
                Initial_Angle(Number_Avalanches) =  estimate_angle(PX(:,t1(na)),PY(:,t1(na)));
                
                Final_Angle(Number_Avalanches)=  estimate_angle(PX(:,t2(na)),PY(:,t2(na)));
                
            end
            
        end
        Avalanche_time{1,nf} = t1;
        Avalanche_time{2,nf} = t2;
        Noavalanches(nf) = Number_Avalanches;
    end
    %fprintf('In file %i the number of avalanches is %i\n',nf,Number_Avalanches);
    
end

%% Resize matrices
Avalanche_duration = Avalanche_duration(1:Number_Avalanches);

Avalanche_particles = Avalanche_particles(1:Number_Avalanches);
Avalanche_displacement = Avalanche_displacement(1:Number_Avalanches);
Avalanche_energy = Avalanche_energy(1:Number_Avalanches);
Avalanche_potential = Avalanche_potential(1:Number_Avalanches);

Normalized_particles = Normalized_particles(:,1:Number_Avalanches);
Normalized_avalanche = Normalized_avalanche(:,1:Number_Avalanches);
Normalized_energy = Normalized_energy(:,1:Number_Avalanches);
Normalized_potential = Normalized_potential(:,1:Number_Avalanches);

mat_particles = mat_particles(:,1:Number_Avalanches);
mat_displacement = mat_displacement(:,1:Number_Avalanches);
mat_energy = mat_energy(:,1:Number_Avalanches);
mat_potential = mat_potential(:,1:Number_Avalanches);

%correlations
correlation_particles = correlation_particles(:,1:Number_Avalanches);
correlation_displacement = correlation_displacement(:,1:Number_Avalanches);
correlation_energy = correlation_energy(:,1:Number_Avalanches);
correlation_potential = correlation_potential(:,1:Number_Avalanches);

% spectrum_particles = spectrum_particles(:,1:Number_Avalanches);
% spectrum_displacement = spectrum_displacement(:,1:Number_Avalanches);
% spectrum_energy = spectrum_energy(:,1:Number_Avalanches);
% spectrum_potential = spectrum_potential(:,1:Number_Avalanches);
Nb_boundary = Nb_boundary(:,1:Number_Avalanches);
DELTAR = DELTAR(1:Number_Avalanches);
Dheight = Dheight(1:Number_Avalanches);
NoParticles_moved = NoParticles_moved(1:Number_Avalanches);
Max_particle_dis = Max_particle_dis(1:Number_Avalanches);
Initial_Angle = Initial_Angle(1:Number_Avalanches);
Final_Angle = Final_Angle(1:Number_Avalanches);
Rotation_step=Rotation_step(:,1:Number_Avalanches);


Displacement_File_nb = Displacement_File_nb(1:Number_Avalanches);
Participation = Participation(1:Number_Avalanches);
In_imafile = In_imafile(1:Number_Avalanches);
Fn_imafile = Fn_imafile(1:Number_Avalanches);
in_trackedfile = in_trackedfile(1:Number_Avalanches);
%% Save results to file
file_save =sprintf('/aline%i/rotdrum%i/o%02d/Avalanches_%i.mat',folder,folder,En,En);


save(file_save,'git_version','maxT','Number_Avalanches','Noavalanches','Avalanche_time', ...
    'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
    'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
    'mat_particles','mat_displacement','mat_energy','mat_potential',...
    'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
    'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
    'Initial_Angle','Final_Angle','Rotation_step','Nb_boundary',...
    'Displacement_File_nb', 'Participation', 'In_imafile','Fn_imafile','in_trackedfile'); 

% save(file_save,'git_version','MaxT','Number_Avalanches','Noavalanches','Avalanche_time', ...
%     'Avalanche_particles','Avalanche_displacement','Avalanche_energy','Avalanche_duration','Avalanche_potential',...
%     'Normalized_particles','Normalized_avalanche','Normalized_energy','Normalized_potential',...
%     'mat_particles','mat_displacement','mat_energy','mat_potential',...
%     'correlation_particles', 'correlation_displacement', 'correlation_energy', 'correlation_potential',...
%     'spectrum_particles','spectrum_displacement','spectrum_energy','spectrum_potential',...
%     'DELTAR','Dheight','NoParticles_moved','Max_particle_dis',...
%     'Initial_Angle','Final_Angle','Rotation_step');


