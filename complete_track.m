
%% Main Sekeleton to track particles

%% Define your folder and Experiment number En and information file of the experiment;
folder = 2;
En = 108;
Cutoff = 11.5;      % minimum peak intensity
MinSep = 6.08;      % minimum separation between peaks 5.5
D = 10;
w = 0.8;


maskfile = '/aline1/rotdrum1/o103/mask103.mat';
load(maskfile,'mk');

avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
load(avanofile,'avan','navfile'); % navfile are the files that contained posible avalanches. 
avalanchefiles = zeros(size(navfile));
NAVFILE = navfile;
save(avanofile,'NAVFILE','-append');
% Find Particle centers using the cluster.


for ii = 1:length(navfile)
    ni = navfile(ii);
    j(ii) = batch(sched,'findparticlecenters',1,...
        {En,ni,folder,D,w,Cutoff,MinSep},'Filedependencies',...
        {'clip.m','chiimg.m','findpeaks.m','ipf.m','discriminate.m'});
    
end

%% Wait for cluster to finish finding particle centers.
for ii = length(navfile)+[-32:0]
   wait(j(ii))
end
clear j
sprintf ('Done finding centers')
% Actualize navfile and separate constinuos files in bunches
for ii = 1:length(navfile)
   ni = navfile(ii);
   fna = (sprintf('/aline%i/rotdrum%i/o%02d/no_avalanche_in_this_file%02d_%05d.mat',folder,folder,En,En,ni));
   avalanchefiles(ii) = exist(fna,'file');
end

navfile = navfile(avalanchefiles == 0);

changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);
save(avanofile,'navfile','-append');

%% Launch mytrack to conect particle positions.
for ii = 1:length(initialfileindex)
   j(ii)=batch(sched,'mytrack',1,{folder,En,ii,initialfileindex(ii),finalfileindex(ii),D,mk},'Filedependencies',{'stickfiles.m','adjacent.m','assignmentoptimal.m'});
     
end
%% Wait for cluster to finish tracking jobs
for ii = length(initialfileindex)+(-32:0)
   wait(j(ii))
end
clear j
sprintf ('Done conecting the centers')
% Find how much particles moves
for ii = 1:length(initialfileindex)
   j(ii)=batch(sched,'displacement',1,{folder,En,ii,initialfileindex(ii),finalfileindex(ii),D});
     
end

%% Get tracked files
 nbtracked = zeros(1,length(initialfileindex));
 for tt =1:length(initialfileindex)
fnn =sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,tt);
nbtracked(tt) = exist(fnn,'file');
 end

nbtracked = find(nbtracked);
save(avanofile,'nbtracked','initialfileindex','finalfileindex','-append');

