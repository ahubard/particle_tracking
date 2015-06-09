%function  [npm]=findava(ni,nframe,En)
%% findinrot  find the particles in the rotating drum.
clf
%File to check 
ofn = sprintf('onestep%2d_%05d.mat',En,ni);
fn=sprintf('onestep%02d_%05d.bin',En,ni);
%% File size
fninfo=sprintf('infof%02d00001.mat',En);
load(fninfo);
Ne=351;
%% Create background
backg=sprintf('back319.mat');
    load(backg,'bg','bbg','bbg2');

%% Ideal Particle
D=10;    %diameter of ideal particle
w=.8;  %Width of ideal particle

%% Create image size grid
%[XX YY]=ndgrid(1:info.Ny,1:info.Nx);
intD = ceil(D);
[xx, yy]=ndgrid(1:intD+1,1:intD+1);
circle=(((xx-(intD/2+1)).^2+(yy-(intD/2+1)).^2)<(intD/2*1.001)^2); %Circle of particle radius. 

%% Parameters for Original image 
% xo=140;    %Drum paramenters
% yo=642;
% R=639;
maskfile = sprintf('/aline/rotdrum/mask430.mat');
load(maskfile);

%Actualize background
iml=readrot(fn,info.Nx,info.Ny,Ne,0,0);

  im=readrot(fn,info.Nx,info.Ny,nframe,0,0);
  bg=max(bg,max(im,iml));
  bbg(iml>0)=max(bbg(iml>0),bg(iml>0)-iml(iml>0)); 
  bbg(im>0)=max(bbg(im>0),bg(im>0)-im(im>0));
 
  first = 1:200;
  second = 201:1076;
  third = 1076:info.Nx;
  
 bbg2=bg;
 bbg2((bg(:,first)-bbg(:,first))<30)=bbg((bg(:,first)-bbg(:,first))<30);
 
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
 
 
 save(backg,'bg','bbg','bbg2');
 sim=clip((bg-im)./bbg2,0,1);
 A=isnan(sim);
 sim(A)=0;

Cutoff=12;      % minimum peak intensity
MinSep=7;      % minimum separation between peaks


% setup for ideal particle

ss=2*fix(D/2+4*w/2)-1;         % size of ideal particle image
os=(ss-1)/2;                   % (size-1)/2 of ideal particle image
[xx, yy]=ndgrid(-os:os,-os:os);  % ideal particle image grid
rr=abs(xx+1i*yy);    % radial coordinate


%Chi image ipf is the ideal image
[ichi]=chiimg(sim,ipf(rr,D,w),[],[],'same');
% find pixel accurate centers using chi-squared

[Np, py, px]=findpeaks(mk./ichi,mk,Cutoff,MinSep);  % find maxima

    simo=sim;
 




sim=clip((bg-iml)./bbg2,0,1);
A=isnan(sim);
sim(A)=0;

pr=(px-xo).^2+(py-yo).^2;   %particle distance to center of drum

pxM=(px>(intD/2)+1);
pxm=(px<(info.Nx-(intD/2)));
pr=(pr<(R-1)^2);
pyM=(py>((intD/2)+1));
pym=(py<(info.Ny-(intD/2)));
ii=find((pxM.*pxm.*pr.*pyM.*pym)>0);
px=px(ii);
py=py(ii);
% simage(sim);
% hold on
% plot(px,py,'.w');
% drawnow;
% hold off



Np=length(ii);
difp=zeros(1,Np);
for np=1:Np
    itx=px(np)-intD/2:px(np)+intD/2;
    ity=py(np)-intD/2:py(np)+intD/2;
     simp=sim(ity,itx);
     simop=simo(ity,itx);
    difim=(simp/mean(simp(:))-simop/mean(simop(:))).*circle;

    difp(np)=sum(abs(difim(:)));
    %difp(np,nframe)=sum(sum(abs((simo-sim).*((((XX-pxs(np,nframe)).^2+(YY-pys(np,nframe)).^2)<5.01^2)))));   %Particle radius is 5
end

npm=sum(difp>30);
 avedif=mean(difp);
 stddif=std(difp);
 
save(ofn,'difp','sim','simo','px','py','avedif','stddif','Np');

