function  D=findparticlecenters(En,ni,folder)
%% findinrot  find the particles in the rotating drum.

fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
load(fn,'IMA');
nfn=sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,ni);

%% File size
fninfo=sprintf('/aline%i/rotdrum%i/o%02d/infof%i00001.mat',folder,folder,En,En);

load(fninfo);
Ne=351;
%% Particle centers
pxs=zeros(3500,Ne);
pys=zeros(3500,Ne);
Npf=zeros(1,Ne);

%% Create background
backg=sprintf('/aline%i/rotdrum%i/back319.mat',folder,folder);
load(backg,'bg','bbg','bbg2');

%% Ideal Particle
D=10;    %diameter of ideal particle
w=0.8;  %Width of ideal particle

%% Create image size grid
%[xx, yy]=ndgrid(1:D+1,1:D+1);
%circle=(((xx-(D/2+1)).^2+(yy-(D/2+1)).^2)<(D/2*1.001)^2); %Circle of particle radius.

%% Parameters for Original image and ideal one
maskfile = sprintf('/aline%i/rotdrum%i/mask430.mat',folder,folder);
load(maskfile);
Cutoff=11.5;      % minimum peak intensity
MinSep=6.08;      % minimum separation between peaks 5.5


% setup for ideal particle

ss=2*fix(D/2+4*w/2)-1;         % size of ideal particle image
os=(ss-1)/2;                   % (size-1)/2 of ideal p
[xx, yy]=ndgrid(-os:os,-os:os);  % ideal particle image grid
rr=abs(xx+1i*yy);    % radial coordinate
%% Main loop over images
for nframe=1:Ne
    
    
    im=IMA(:,:,nframe);
    bg=max(bg,im);
    bbg(im>0)=max(bbg(im>0),bg(im>0)-im(im>0));
    
    first = 1:200;
    second = 201:1076;
    third = 1076:info.Nx;
    
    bbg2=bg;
    bbg2((bg(:,first)-bbg(:,first))<40)=bbg((bg(:,first)-bbg(:,first))<40);
    
    bgp = bg(:,second);
    bbgp = bbg(:,second);
    bbg2p = bbg2(:,second);
    bbg2p((bgp-bbgp)<50)= bbgp((bgp-bbgp)<50);
    bbg2(:,second) = bbg2p;
    
    bgp = bg(:,third);
    bbgp = bbg(:,third);
    bbg2p = bbg2(:,third);
    bbg2p((bgp-bbgp)<25)= bbgp((bgp-bbgp)<25);
    bbg2(:,third) = bbg2p;
    
    
    
    low = 0.1;
    high = 0.98;
    sim=(clip((bg-im)./bbg2,low,high)-low)/(high-low);
    A=isnan(sim);
    sim(A)=0;
    

    
    %Chi image ipf is the ideal image
    [ichi]=chiimg(sim,ipf(rr,D,w),[],[],'same');
    % find pixel accurate centers using chi-squared
    
    [Np, py, px]=findpeaks(mk./ichi,mk,Cutoff,MinSep);  % find maxima
    
    pys(1:Np,nframe)=py;
    pxs(1:Np,nframe)=px;
    Npf(nframe)=Np;
end
%save(backg,'bg','bbg','bbg2');
save(nfn,'pys','pxs','Npf');


end


