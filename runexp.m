%% runexp run de experiment, finds the files where there are avalanches and
%then create .mat file with the images, the particles positions between the
%first and last image. And finally erases the bin files to keep only .mat
%where there was some avalanche. 

%% Runs experiment
En=15;
ifwrong=1;
while(ifwrong)
avasnap
fprintf('Done taking data');

%% Check there was no error
ifwrong=(max(checking(1:ns-1)-((1:ns-1)*NS)-NAC1)>info.nbuf);
end
%% loads info file and searchs for avalanches that last more than one file
Avanostep=sprintf('Avanonestep%02d.mat',En);
load(Avanostep);
dnr=diff(avan(2,1:na));    %Difference between rotation in files
ia=find(dnr==0)+1;           %If dif is zero then there was an avalanche
navfile=unique([ia-1 ia]);      %Add next file where the avalanche is still ocurring
%% Find small avalanches that happened in one single file.
mm=zeros(1,na);
mm(navfile)=1;
serchn=find(mm==0);  %Files with posibles short duration(one file) avalan
SE=length(serchn);
nframe=1;
NPM=zeros(1,SE);
for n=1:SE
    ni=serchn(n);
    [npm]=findava(ni,nframe,En);  %Number of particles per frame that moved more than rotation
    NPM(n)=npm;
end
hh=find(NPM>0);     %Files with short avalanches ocurred. 
navfile=sort([navfile hh]);
save(Avanostep,'navfile','-append');
%% Go over navfile bin files and create mat files while tracking for the centers positions. 
fprintf('Writing .mat files');


for n=1:length(navfile);
ni=navfile(n);
fnb=sprintf('onestep%02d%05d.bin',En,ni);
fn=sprintf('onestep%02d%05d.mat',En,ni);
IMA=zeros(info.Ny,info.Nx,NS+1);
for nframe=1:NS+1
im=readrot(fnb,info.Nx,info.Ny,nframe,0,0);
IMA(:,:,nframe)=im;
end
save(fn,'IMA');

end
 %% Erase bin files
% fprintf('Erasing .bin files');
% for ni=1:na
%  fnb=sprintf('onestep13%05d.bin',ni);
%  delete(fnb);
% end
%% Track particles on navfiles
fprintf('Tracking Particles')
%load background
backg=sprintf('back13.mat');
    load(backg,'bg','bbg','bbg2');
% Ideal Particle
D=10;    %diameter of ideal particle
w=.9;  %Width of ideal particle

% Create image size grid
[xx, yy]=ndgrid(1:D+1,1:D+1);
circle=(((xx-(D/2+1)).^2+(yy-(D/2+1)).^2)<(D/2*1.001)^2); %Circle of particle radius. 

% Parameters for Original image and ideal one
load('mask13.mat');
Cutoff=9;      % minimum peak intensity
MinSep=6;      % minimum separation between peaks


% setup for ideal particle

ss=2*fix(D/2+4*w/2)-1;         % size of ideal particle image
os=(ss-1)/2;                   % (size-1)/2 of ideal particle image
[xx, yy]=ndgrid(-os:os,-os:os);  % ideal particle image grid
rr=abs(xx+1i*yy);    % radial coordinate
   
for n=1:length(navfile);
ni=navfile(n);
findparticlecenters;
end
    