%function [keep] = discriminate(folder,En,ni,D,w,Cutoff,MinSep)

function [keep, IMA,mk,bk1,bk2,info,maxdifp,participationratio,Num_p, standardev] = discriminate(folder,En,ni,D,w,Cutoff,MinSep)
%% Reads file in folder"aline<folder>/rotdrum<folder>/o<En>/onestep<En>_<ni>...
%and checks if there is an avalanche in it. Using that D is the diameter of
%ideal particle,and w the ideal particle width of the brightness. 
% keep = 0, no avalanche or perceptible movement, keep = 1, only one of the
% conditions was fullfilled, probably there is just solid angle rotation,
% keep = 2, both conditions fullfilled, some of the particles moved in an
% avalanche. IMA is the series of pictures with images, bk1,bk2 are the
% background, mk the mask. 
% Cutoff=11.5;      % minimum peak intensity
% MinSep=6.08;      % minimum separation between peaks 5.5
%% Cutoff to decide if there is an avalanche. 
cutoffpr = 8e-4;
cutoffdiff = 10;
cutoffstd = 6;
%% Check if files exist and load them 
if (En >100)  
     bgfile = sprintf('/aline2/rotdrum2/o105/back105.mat');
else
     bgfile = sprintf('/aline%i/rotdrum%i/back_%ito%i.mat',folder,folder,15,19);
end
load(bgfile,'bk1','bk2');

bkk1 = bk1;
bkk2 = bk2;



avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
gotmaskinfo = whos(matfile(avanofile),'mk','xo','yo','R');

bgfile = sprintf('/aline%i/rotdrum%i/o%i/back%i.mat',folder,folder,En,En);

if (En > 100)
    fno = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
    
    if (size(gotmaskinfo,1) < 4)
%         comment1 = char('used /aline1/rotdrum1/o103/mask103.mat as a mask');
%         save(fno,'comment1','-append');
        avanofile = '/aline1/rotdrum1/o103/mask103.mat';
    end
    
    if (exist(bgfile,'file') == 0)
%         comment2 = char('used /aline2/rotdrum2/o105/back105.mat as a background');
%         save(fno,'comment2','-append');
        bgfile = '/aline2/rotdrum2/o105/back105.mat';
    end
else
    fno=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d%05d.mat',folder,folder,En,En,ni);
    
    if (size(gotmaskinfo,1) < 4)
%         comment1 = char('used /aline1/rotdrum1/mask13.mat as a mask');
%         save(fno,'comment1','-append');
        avanofile = '/aline1/rotdrum1/mask13.mat';
    end
    
    if (exist(bgfile,'file') == 0)
%         comment2 = char('used /aline2/rotdrum2/back20.mat as a background');
%         save(fno,'comment2','-append');
        bgfile = '/aline2/rotdrum2/back20.mat';
    end
end



load(avanofile,'mk','xo','yo','R');

load(bgfile,'bk1','bk2');
bk1 = max(bk1,bkk1);
bk2 = max(bk2,bkk2);


bg = bk1*0;
bbg = bg;
bbg2 =bbg;
load(fno,'IMA');


info = struct('Nx',size(IMA,2),'Ny',size(IMA,1),'Numframe',size(IMA,3));

%% Find particles positions of first image

% Create filter
se = strel('disk',D+4);
filtercof = 50;

% setup for ideal particle

ss = 2*fix(D/2+4*w/2)-1;               % size of ideal particle image
os = (ss-1)/2;                         % (size-1)/2 of ideal p
[xx, yy] = ndgrid(-os:os,-os:os);      % ideal particle image grid
rr = abs(xx+1i*yy);                    % radial coordinate

%Normalize first image by amount of light
im = IMA(:,:,1);
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

simo = sim; 
%Chi image ipf is the ideal image
[ichi] = chiimg(sim,ipf(rr,D,w),[],[],'same');
% find pixel accurate centers using chi-squared
[Num_p, py, px] = findpeaks(mk./ichi,mk,Cutoff,MinSep);  % find maxima
%% Compare first image with last image at particles positions. 

%Normalize last image
iml = IMA(:,:,info.Numframe);
bgmask = (imopen(iml.*mk./bk1,se))>filtercof;
normcoef = sum(sum(iml.*bgmask))/sum(bgmask(:));
iml = iml/normcoef;

sim = clip((bg-iml)./bbg2,0,1);
A = isnan(sim);
sim(A) = 0;


%Create mask of circles to put around the particles
intr = ceil(D/2);
[xx, yy] = ndgrid(1:2*intr+1,1:2*intr+1);
circle = (((xx-(intr+1)).^2+(yy-(intr+1)).^2)<(intr*1.001)^2); %Circle of particle radius.



pr = (px-xo).^2+(py-yo).^2;   %particle distance to center of drum
pxM = (px > (intr+1));          %Get rid of the edges.
pxm = (px < (info.Nx-intr));
pr = (pr < (R-1)^2);
pyM = (py > (intr+1));
pym = (py < (info.Ny-intr));
ii = find((pxM.*pxm.*pr.*pyM.*pym)>0);
px = px(ii);
py = py(ii);



Np=length(ii);
difp=zeros(1,Np);


for np=1:Np
    
    itx = px(np)-intr:px(np)+intr;
    ity = py(np)-intr:py(np)+intr;
    simp = sim(ity,itx);
    simop = simo(ity,itx);
    difim = (simp/mean(simp(:))-simop/mean(simop(:))).*circle;
    difp(np) = sum(abs(difim(:)));
    
end

maxdifp = max(difp);
participationratio =  sum(difp.^4)/(sum(difp.^2)).^2;
standardev = std(reshape(sim-simo,1,[]));
keep = ((maxdifp > cutoffdiff) + (participationratio > cutoffpr) + ...
    (standardev < cutoffstd) );


