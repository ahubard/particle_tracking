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
surfacebestfitline = zeros(2,Ne);
%% Create background
bgfile = sprintf('/aline%i/rotdrum%i/o%i/back%i.mat',folder,folder,En,En);
load(bgfile,'bk1','bk2');
bg = bk1*0;
bbg = bg;
bbg2 =bbg;
%% Ideal Particle
D=10;    %diameter of ideal particle
w=0.8;  %Width of ideal particle

%% Create filter
se = strel('disk',D+4);
filtercof = 50;

%% Parameters for Original image and ideal one
maskfile = sprintf('/aline%i/rotdrum%i/mask%i.mat',folder,folder,En,En);
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
      
    
    
      
    
    im = IMA(:,:,nframe);
    bgmask = (imopen(im.*mk./bk1,se))>filtercof;
    normcoef = sum(sum(im.*bgmask))/sum(bgmask(:));
    
    im = im/normcoef;
    
  
    bg=max(bk1,im);
    bbg(im>0)=max(bk2(im>0),bg(im>0)-im(im>0));
    bbg2(bgmask) = bg(bgmask); 
    bbg2(bgmask==0) = bbg(bgmask==0);
    
    
    low = 0.1;
    high = 0.95;
    sim=(clip((bg-im)./bbg2,low,high)-low)/(high-low);
    A=isnan(sim);
    sim(A)=0;
    
    
    %Chi image ipf is the ideal image
    [ichi]=chiimg(sim,ipf(rr,D,w),[],[],'same');
    % find pixel accurate centers using chi-squared
    
    [~, py, px]=findpeaks(mk./ichi,mk,Cutoff,MinSep);  % find maxima
    
% Keep only insiders 
    [xs ix] = sort(ceil(px/binsize));
    bindex = [0 ;find(diff(xs)>0)];
    ys = py(ix);
    ya = zeros(1,length(bindex)-1);
    xa = ya;
    
    for bin = 1:length(bindex)-1;
        ya(bin) = min(ys(bindex(bin)+1:bindex(bin+1)));
        xa(bin) = binsize*(bin + xs(1)-1);
    end
    
   linearfit = polyfit(xa,ya,1);
   nvector = [-linearfit(1) 1]/(sqrt(linearfit(1)^2+1)); %normalvector to bestfit line
   lo = linearfit(2)*nvector(2); %Distance of bestfitline to origin. 
   dpointtoline = px*nvector(1) + py*nvector(2) - lo;
   insiders = (dpointtoline > -D*4); %Points that are in the region of interest
   surfacebestfitline(:,nframe) = linearfit;
   
   %Get rid of outliers
   Npf(nframe) = sum(insiders);
   pxs(1:Npf(nframe),nframe) = px(insiders);
   pys(1:Npf(nframe),nframe) = py(insiders);
   
   
    
    
   
end
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
save(nfn,'git_version','pys','pxs','Npf','surfacebestfitline');


end


