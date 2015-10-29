function improve_background(folder,En)

sprintf('If you wanna to get the background of folder %i and experiment series %i press enter, otherwise cancel',folder,En)
pause()
sprintf('OK, here we go, we will start by opening the scheduler.')
%% Main Sekeleton to track particles
sched=findResource('scheduler', 'type', 'jobmanager', 'lookupURL','poincare.engr.ccny.cuny.edu');  %Open scheduler

[git_version, ~] = evalc('system(''git describe --dirty --alway'')')


%% Experiment file
avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
load(avanofile,'avan','navfile','mk'); % navfile are the files that contained posible avalanches.

%% Get background
Nbfiles = length(navfile);
bgfile = sprintf('/aline%i/rotdrum%i/o%i/back%i.mat',folder,folder,En,En);

schedulerindex = 1;
for ii = 1:Nbfiles
    ni = navfile(ii);
    jj(schedulerindex) = batch(sched,'getbackground',1,{En,ni,folder,1});
    schedulerindex = schedulerindex+1;
end

for ii = (-32:-1)+schedulerindex
    wait(jj(ii));
end

bk = mk*0;

for ii = 1:Nbfiles
    ni = navfile(ii);
    
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
for ii = 1:Nbfiles
    ni = navfile(ii);
    jj(schedulerindex) = batch(sched,'getbackground',1,{En,ni,folder,2});
    schedulerindex = schedulerindex+1;
end

for ii = (-32:-1)+schedulerindex
    wait(jj(ii));
end

bk = mk*0;

for ii = 1:Nbfiles
    ni = navfile(ii);
    if (En > 100)
        fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
    else
        fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d%05d.mat',folder,folder,En,En,ni);
    end
    
    load(fn,'bk2');
    bk = max(bk,bk2);
end

bk2 = bk;
save (bgfile,'bk2','git_version','-append');
sprintf('Done with second part of the background')
