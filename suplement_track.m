%function suplement_track(folder,En)

sprintf('If you wanna complete the tracking of  folder %i and experiment series %i again press enter, otherwise cancel',folder,En)
pause()
sprintf('OK, here we go, we will start by opening the scheduler.')
%% Main Sekeleton to track particles
sched=findResource('scheduler', 'type', 'jobmanager', 'lookupURL','poincare.engr.ccny.cuny.edu');  %Open scheduler

[git_version, ~] = evalc('system(''git describe --dirty --alway'')')

%% Define your folder and Experiment number En and information file of the experiment;
% folder = 2;
% En = 109;


D = 10;
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);


%% Experiment file
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile,'avan','navfile','mk'); % navfile are the files that contained posible avalanches.

changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);


schedulerindex =1; 
try_again = zeros(1,length(initialfileindex));
for ii = 1:length(initialfileindex)
    filekernel =sprintf('Tracked_%i',ii);
    fna = (sprintf('%sWarning_The_file_%s_does_not_exist.mat',filedirectory,filekernel));
    try_again(ii) = exist(fna,'file');
end
numtrack = find(try_again);
sprintf('I need to reconect %i files',length(numtrack))
%% Launch mytrack to conect particle positions.
%for ij = 1:length(numtrack)
for ij = 1:10
    ii = numtrack(ij);
   jj(schedulerindex)=batch(sched,'mytrack',1,{folder,En,ii,initialfileindex(ii),finalfileindex(ii),D,mk},'Filedependencies',{'stickfiles.m','adjacent.m','assignmentoptimal.m'});
   schedulerindex = schedulerindex+1; 
end
%% Wait for cluster to finish tracking jobs
% for ii = (-32:-1)+schedulerindex
%     wait(jj(ii));
% end

%sprintf ('Done conecting the centers,do displacement function now.')

% Find how much particles moves
%count = displacementscript(folder,En);

% %% Get tracked files
%  nbtracked = zeros(1,length(initialfileindex));
%  for tt =1:length(initialfileindex)
% fnn =sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,tt);
% nbtracked(tt) = exist(fnn,'file');
%  end
% 
% nbtracked = find(nbtracked);
% save(avanofile,'nbtracked','initialfileindex','finalfileindex','-append');

