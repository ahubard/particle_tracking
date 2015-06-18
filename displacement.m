% Run function over the NEn = NItracked(ii) files that are already tracked. 

%Get tracked files;
% yestrack = zeros(1,length(initialfileindex)); initialfileindex comes from
% stickfiles
% for nbfile = 1 :length(initialfileindex)
%     TRACKFILE =sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,nbfile);
%     yestrack(nbfile) = exist(TRACKFILE,'file');
% end
%  
% NItracked = find(yestrack);

%!git add displacement.m
%!git commit -m "added git version in output"



function [participationratio]=displacement(folder,En,NEn,initial,final,D)

%Finds which particles move more than threshold and their diference over
%time

%% File to load
%infofile = sprintf('/aline%i/rotdrum%i/o%02d/infof%02d%05d.mat',folder,folder,En,En,1);
maskfile = sprintf('/aline%i/rotdrum%i/mask430.mat',folder,folder);
backg = sprintf('/aline%i/rotdrum%i/back319.mat',folder,folder);
fni = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,initial);
fnf = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,final);
fnt =sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,NEn);
load(fni,'IMA');
imi = IMA(:,:,1);
load(fnf, 'IMA');
imf = IMA(:,:,351);
load(fnt,'PX','PY');
numframes = size(PX,2);

%% Load background and general mask

load(backg,'bg','bbg','bbg2');
load(maskfile);
%% IMAGES info

xsize=size(IMA,2);
ysize=size(IMA,1);


%% Get the particles involved in avalanche.
iframe1=(PX(:,1)>0);                             %Particles indexes on first frame
ilastframe=(PX(:,numframes)>0);                         %On last frame
iframe=find((ilastframe.*iframe1)==1);          %Particles that didnt dissapear in the tracking
nbparticles=length(iframe);

px1=PX(iframe,1);
py1=PY(iframe,1);
%%  Get image difference at particles positions between first and last image

%Normalized first image

firstimage=clip((bg-imi)./bbg2,0,1);
A=isnan(firstimage);
firstimage(A)=0;

%Normalized last image
lastimage=clip((bg-imf)./bbg2,0,1);
A=isnan(lastimage);
lastimage(A)=0;

imagediff = mk.*(lastimage-firstimage).^2;
centerfraction = sum(sum(imagediff(100:300,400:800)))/sum(imagediff(:));
%simage(imagediff);
%drawnow;
%% Check in which particle locations image changed
IDIFF=zeros(1,nbparticles);
ISTD=zeros(1,nbparticles);
cutoff=0.0015;
for nf=1:nbparticles
    nfxinitial=max(px1(nf)-D/2,1);
    nfxfinal=min(px1(nf)+D/2,xsize);
    nfyinitial=max(py1(nf)-D/2,1);
    nfyfinal=min(py1(nf)+D/2,ysize);
    
    nfdiff=imagediff(nfyinitial:nfyfinal,nfxinitial:nfxfinal);
    IDIFF(nf)=mean(nfdiff(:));
    ISTD(nf)=std(nfdiff(:));
end

   diskmove=iframe(IDIFF>cutoff);
  
   






%% Find displacements of such particles btw nframes
increment=10;
%dt = increment/2;
dx =(PX(diskmove,1+increment:1:numframes)-PX(diskmove,1:1:numframes-increment))/increment;
dy =(PY(diskmove,1+increment:1:numframes)-PY(diskmove,1:1:numframes-increment))/increment;
dr = dx.^2+dy.^2;
%dr=sqrt(dr);
totaltr=((PX(diskmove,numframes)-PX(diskmove,1)).^2+(PY(diskmove,numframes)-PY(diskmove,1)).^2);

%totdisplacement=sum(totaltr);
participationratio=sum(totaltr.^4)/(sum(totaltr.^2)^2);
%% Save results
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
fnn =sprintf('/aline%i/rotdrum%i/o%02d/Trackeds_%i.mat',folder,folder,En,NEn);
save(fnn,'git_version','centerfraction','dr','diskmove','increment','participationratio','imagediff','folder','En','NEn','initial','final','D');
