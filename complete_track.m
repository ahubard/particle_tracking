function complete_track(folder,En)

sprintf('If you wanna track folder %i and experiment series %i press enter, otherwise cancel',folder,En)
pause()
sprintf('OK, here we go, we will start by opening the scheduler.')
%% Main Sekeleton to track particles
sched=findResource('scheduler', 'type', 'jobmanager', 'lookupURL','poincare.engr.ccny.cuny.edu');  %Open scheduler

[git_version, ~] = evalc('system(''git describe --dirty --alway'')')
%% Check files with images 
files_index = review_files(folder,En);
avalanchefiles = zeros(size(files_index));
%% Define your folder and Experiment number En and information file of the experiment;

Cutoff = 11.5;      % minimum peak intensity
MinSep = 6.08;      % minimum separation between peaks 5.5
D = 10;
w = 0.8;

%% Experiment file
avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
%load(avanofile,'avan','mk'); 
load(avanofile); 
schedulerindex = 1;

%% Get background
Nbfiles = length(files_index);
bgfile = sprintf('/aline%i/rotdrum%i/o%i/back%i.mat',folder,folder,En,En);
jumps = ceil(Nbfiles/64);

bgfileexists = whos(matfile(bgfile));

if (length(bgfileexists) < 1)
    
    for ii = 1:jumps:Nbfiles
        ni = files_index(ii);
        jj(schedulerindex) = batch(sched,'getbackground',1,{En,ni,folder,1});
        schedulerindex = schedulerindex+1;
    end
    
    for ii = (-32:-1)+schedulerindex
        wait(jj(ii));
    end
    
    bk = mk*0;
    
    for ii = 1:jumps:Nbfiles
        ni = files_index(ii);
        
        if (En > 100)
            fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
        else
            fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d%05d.mat',folder,folder,En,En,ni);
        end
        
        load(fn,'bk1');
        bk = max(bk,bk1);
    end
    
    bk1 = bk;
    save (bgfile,'bk1');
    
    sprintf('Done with first part of the background')
    
    if(length(bgfileexists) < 2)
        for ii = 1:jumps:Nbfiles
            ni = files_index(ii);
            jj(schedulerindex) = batch(sched,'getbackground',1,{En,ni,folder,2});
            schedulerindex = schedulerindex+1;
        end
        
        for ii = (-32:-1)+schedulerindex
            wait(jj(ii));
        end
        
        bk = mk*0;
        
        for ii = 1:jumps:Nbfiles
            ni = files_index(ii);
            if (En > 100)
                fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
            else
                fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d%05d.mat',folder,folder,En,En,ni);
            end
            
            load(fn,'bk2');
            bk = max(bk,bk2);
        end
        
        bk2 = bk;
        save (bgfile,'bk2','-append');
        sprintf('Done with second part of the background')
    end
end


%% Find Particle centers using the cluster.
%schedulerindex = 2*length(1:jumps:Nbfiles)+1;


for ii = 1:Nbfiles
    ni = files_index(ii);
    jj(schedulerindex) = batch(sched,'findparticlecenters',1,...
        {En,ni,folder,D,w,Cutoff,MinSep},'Filedependencies',...
        {'clip.m','chiimg.m','findpeaks.m','ipf.m','discriminate.m'});
    schedulerindex = schedulerindex+1;
    
end
sprintf('Be patient, Im trying to find the centers') 

%% Wait for cluster to finish finding particle centers.
for ii = (-32:-1)+schedulerindex
    wait(jj(ii));
end

sprintf ('Done finding centers')

% Actualize navfile and separate continuos files in bunches
for ii = 1:Nbfiles
    ni = files_index(ii);
    fna = (sprintf('/aline%i/rotdrum%i/o%02d/no_avalanche_in_this_file%02d_%05d.mat',folder,folder,En,En,ni));
    avalanchefiles(ii) = exist(fna,'file');
end

navfile = files_index(avalanchefiles == 0);

changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);
 save(avanofile,'navfile','-append');

sprintf('There is %i files that contain avalanches',length(navfile))

%% Launch mytrack to conect particle positions.
for ii = 1:length(initialfileindex)
   jj(schedulerindex)=batch(sched,'mytrack',1,{folder,En,ii,initialfileindex(ii),finalfileindex(ii),D,mk},'Filedependencies',{'stickfiles.m','adjacent.m','assignmentoptimal.m'});
   schedulerindex = schedulerindex+1; 
end

%% Wait for cluster to finish tracking jobs
for ii = (-32:-1)+schedulerindex
    wait(jj(ii));
end

sprintf ('Done conecting the centers')
% Find how much particles moves
count = displacementscript(folder,En);

Number_Avalanches = avalanche_size(folder,En);
sprintf('For experiment %i I found %i avalanches in %i different files',En,Number_Avalanches,count)

% %% Get tracked files
%  nbtracked = zeros(1,length(initialfileindex));
%  for tt =1:length(initialfileindex)
% fnn =sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,tt);
% nbtracked(tt) = exist(fnn,'file');
%  end
% 
% nbtracked = find(nbtracked);
% save(avanofile,'nbtracked','initialfileindex','finalfileindex','-append');

