
%% Main Sekeleton to track particles

%% Define your folder and Experiment number En and information file of the experiment;
folder = 2;
En = 105;
Cutoff = 11.5;      % minimum peak intensity
MinSep = 6.08;      % minimum separation between peaks 5.5
D = 10;
w = 0.8;

avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
load(avanofile,'avan','navfile'); % navfile are the files that contained posible avalanches. 
avalanchefiles = zeros(size(navfile));
NAVFILE = navfile;
save(avanofile,'-append');
%% Find Particle centers using the cluster.

for ii = 1:length(navfile)
    ni = navfile(ii);
    jobs(ii) = batch(sched,'findparticlecenters',1,...
        {En,ni,folder,D,w,Cutoff,MinSep},'Filedependencies',...
        {'clip.m','chiimg.m','findpeaks.m','ipf.m','discriminate.m'});
    
end

%% Wait for cluster to finish finding particle centers.
for ii = length(navfile)+[-32:0]
    wait(jobs(ii))
end
clear jobs
sprintf ('Done finding centers')
%% Actualize navfile and separate constinuos files in bunches
for ii = 1:length(navfile)
    ni = navfile(ii);
    fna = (sprintf('/aline%i/rotdrum%i/o%02d/no_avalanche_in_this_file%02d_%05d.mat',folder,folder,En,En,ni));
    avalanchefiles(ii) = exist(fna,'file');
end

navfile = navfile(avalanchefiles == 0);

changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);

%% Launch mytrack to conect particle positions.
for ii = 1:length(initialfileindex)
    trjob(ii)=batch(sched,'mytrack',1,{folder,En,ii,firstfile(ii),finalfile(ii),D},'Filedependencies',{'stickfiles.m','adjacent.m','assignmentoptimal.m'});
     
end
%% Wait for cluster to finish tracking jobs
for ii = length(navfile)+(-32:0)
    wait(trjobs(ii))
end
clear trjobs
sprintf ('Done conecting the centers')
%% Find how much particles moves
for ii = 1:length(initialfileindex)
    trjob(ii)=batch(sched,'displacement',1,{folder,En,ii,initialfileindex(ii),finalfileindex(ii),D});
     
end

