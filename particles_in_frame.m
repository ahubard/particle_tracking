function [px, py, Nbp, sim] = particles_in_frame(im,filedirectory,En,D,w)

Cutoff = 11.5;      % minimum peak intensity
MinSep = 6.08;      % minimum separation between peaks 5.5





%load background info
bgfile = sprintf('%sback%i.mat',filedirectory,En);
load(bgfile, 'bk1','bk2');

%load mask info
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile,'mk','xo','yo','R');

bg = bk1*0;
bbg = bg;
bbg2 =bbg;

% Create filter
se = strel('disk',D+4);
filtercof = 50;

% setup for ideal particle

ss = 2*fix(D/2+4*w/2)-1;               % size of ideal particle image
os = (ss-1)/2;                         % (size-1)/2 of ideal p
[xx, yy] = ndgrid(-os:os,-os:os);      % ideal particle image grid
rr = abs(xx+1i*yy);                    % radial coordinate

%Normalize first image by amount of light
bgmask = (imopen(im.*mk./bk1,se))>filtercof;
normcoef = sum(sum(im.*bgmask))/sum(bgmask(:));

im = im/normcoef;

%Create background
bg = max(bk1,im);
bbg(im>0) = max(bk2(im>0),bg(im>0)-im(im>0));
bbg2(bgmask) = bg(bgmask);
bbg2(bgmask==0) = bbg(bgmask==0);

%Substract background
low = 0.1;
high = 0.95;
sim = (clip((bg-im)./bbg2,low,high)-low)/(high-low);
A =isnan(sim);
sim(A) = 0;

%Chi image ipf is the ideal image
[ichi] = chiimg(sim,ipf(rr,D,w),[],[],'same');
% find pixel accurate centers using chi-squared
[Nbp, py, px] = findpeaks(mk./ichi,mk,Cutoff,MinSep);  % find maxima
sim = sim.*mk;