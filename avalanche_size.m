function Number_Avalanches = avalanche_size(folder,En)
%% Gather the basic info

avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
%avanofile = sprintf('Avanonestep%i.mat',En);

if (exist(avanofile,'file'))
    load (avanofile,'count');
else
    save(sprintf('Warning. The file: %s does not exist.mat',avanofile));
    error('Error, avanofile does not exist');
end

[git_version, ~] = evalc('system(''git describe --dirty --alway'')');

%% Initialize variables.
Avalanche_size = zeros(1,4*count);
Avalanche_duration = zeros(1,4*count);
Normalized_avalanche = zeros(101,4*count);
Number_Avalanches = 0;
%% Cutoff values
cutoffperparticle = 0.01;    %If drfil (filtered particle displacement) is smaller
% cutof then the particle didn move.
cutoffofsum = 0.01;
windowSize = 5;
minstepsbtwav = 10; %same as increment. Minimal amount of steps whit no displacement to say the avalanche is over.
b = (1/minstepsbtwav)*ones(1,minstepsbtwav);
a = 1;


%% Main loop
for nf = 1:count
    fnn =sprintf('/aline%i/rotdrum%i/o%02d/Displacement_%i.mat',folder,folder,En,nf);
    %fnn =sprintf('Displacement_%i.mat',nf);
    clear('diskmove','drraw','drfil');
    
    if(exist(fnn,'file'))
        load(fnn,'diskmove','drfil','drraw');
    else
        save(sprintf('Warning. The file: %s does not exist.mat',fnn));
        error('Error, Displacementfile does not exist');
    end
    
    if (diskmove)%Check there was no error in the file. 
        
        drmask = (drfil>=cutoffperparticle);  %Use filter dr to determine which particles moved at each time step.
        avalanche = sum(drraw(:,ceil(windowSize/2):size(drraw,2)-floor(windowSize/2)).*drmask);
        findavalanche = filter(b,a,(avalanche>cutoffofsum)); %use both numbers tdecide, if they are moving they should move for at least 10 time steps, if they stop the same
        % Find how many avalanches per file findavalanche = 0 means no
        % avalanche there.
        
        %Posible beginings and endings of avalanches
        t1 = find(diff(findavalanche>0)==1);  %Time steps avalanche begins
        if(isempty(t1)) 
            t1 = 1;
        end
        
        t2 = min(find(diff(findavalanche>0)==-1)+1,length(avalanche)); %Time steps avalanche ends.
        if(isempty(t2))
            t2 = length(avalanche);
        end
        if (t1(1)>= t2(1)) %In case avalanche starts at first step.
            t1 = [1 t1];
        end
        
        if(t2(end)<=t1(end))%In case avalanche finish at the last frame of the file.
            t2 = [t2 length(avalanche)];
        end
        
        
        for na = 1:length(t1)
            if( sum((findavalanche(t1(na):min(t2(na),length(avalanche))))>=1-eps)) %Check if there is indeed avalanches between t1-t2
                Number_Avalanches = Number_Avalanches+1;
                Avalanche_size(Number_Avalanches) = sum (avalanche(t1(na):t2(na)));
                Avalanche_duration(Number_Avalanches) = length(t1(na):t2(na));
                avalanchenormalized = interp1( (0:(t2(na)-t1(na)))/(t2(na)-t1(na)),avalanche(t1(na):t2(na)),(0:.01:1),'pchip');
                Normalized_avalanche(:,Number_Avalanches) = avalanchenormalized/max(avalanchenormalized);
            end
            
        end
    end
end
Avalanche_size = Avalanche_size(1:Number_Avalanches);
Avalanche_duration = Avalanche_duration(1:Number_Avalanches);
Normalized_avalanche = Normalized_avalanche(:,1:Number_Avalanches);

file_save =sprintf('/aline%i/rotdrum%i/o%02d/Avalanches_%i.mat',folder,folder,En,En);
%file_save =sprintf('Avalanches_%i.mat',En);
save(file_save,'git_version','Number_Avalanches','Avalanche_duration','Avalanche_size','Normalized_avalanche');
    
    