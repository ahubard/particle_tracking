function [px,py,NPF,initial] = stickfiles(folder,En,initial,final,D)


fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,initial);



load(fn,'pxs','pys','Npf');


nfiles = final-initial+1;

    px = zeros(3500,nfiles*350);
    py = zeros(3500,nfiles*350);
    NPF = zeros(1,nfiles*350);
    Npf = 3500*ones(1,351);
    
    
    
    for ii = initial:final
        fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,ii);
        load(fn,'pxs','pys','Npf');
        [i0 t0] = find(pxs == 0);
        [t0 it] = unique(t0,'first');
        Npf(t0) = i0(it)-1;
        MNpf = max(Npf);
        px(1:MNpf,((ii - initial)*350)+1:(ii-initial+1)*350) = pxs(1:MNpf,1:350);
        py(1:MNpf,((ii - initial)*350)+1:(ii-initial+1)*350) = pys(1:MNpf,1:350);
        NPF(((ii - initial)*350)+1:(ii-initial+1)*350) = Npf(1:350);
    end
    

%This is how to get initial and final.

% avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En)
% load(avanofile,'avan','navfile');
% changefileindex = find(diff(avan(1,navfile))>1);
% initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
% finalfileindex = navfile(changefileindex);