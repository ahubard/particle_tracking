function  D = refinebackground(En,ni,folder)
%% refinebackground
% uses tracked particles in previous run to get a new background that does
% not over expose regions. 
%% LOADPOSITIONS AND IMAGE
fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
load(fn,'IMA');
nfn=sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,ni);
load(nfn,'pxs','pys','Npf');
fninfo=sprintf('/aline2/rotdrum2/infof10600001.mat');
load(fninfo);
load(sprintf('/aline%i/rotdrum%i/back319.mat',folder,folder),'BG');


D = 10;
R = ceil(D/2)+1;
Ne=351;
%% background and normalization matrices

backni2 = IMA(:,:,1)*0;
normni2 = backni2;


%% Main loop over images
for nframe=1:Ne

     mask = circlesmask(pxs(1:Npf(nframe),nframe),pys(1:Npf(nframe),nframe),info.Nx,info.Ny,Npf(nframe),R);
     im=IMA(:,:,nframe);
     backni2 = backni2 + mask.*(BG-im);
     normni2 = normni2 + mask;
     
     

end

save(fn,'backni2','normni2','-append');


end
