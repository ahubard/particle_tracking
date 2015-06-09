function [px,py,NPF,initial] = stickfiles(folder,En,initial,final,D)
cutoffpr = 8e-4;
cutoffdiff = 10;

%% Make sure first file moved more than solid angle.
infofile = sprintf('/aline%i/rotdrum%i/o%02d/infof%02d%05d.mat',folder,folder,En,En,1);
maskfile = sprintf('/aline%i/rotdrum%i/mask430.mat',folder,folder);
backg = sprintf('/aline%i/rotdrum%i/back319.mat',folder,folder);
fno = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,initial);
fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,initial);

load(infofile);
load(maskfile);
load(backg,'bg','bbg','bbg2');
load(fno,'IMA');
load(fn,'pxs','pys','Npf');

px = pxs(:,1);
py = pys(:,1);
NPF = length(px);
im = IMA(:,:,1);
iml = IMA(:,:,351);


intD = ceil(D);
[xx, yy]=ndgrid(1:intD+1,1:intD+1);
circle=(((xx-(intD/2+1)).^2+(yy-(intD/2+1)).^2)<(intD/2*1.001)^2); %Circle of particle radius.

simo=clip((bg-im)./bbg2,0,1);
A=isnan(simo);
simo(A)=0;

sim=clip((bg-iml)./bbg2,0,1);
A=isnan(sim);
sim(A)=0;


pr=(px-xo).^2+(py-yo).^2;   %particle distance to center of drum

pxM=(px>(intD/2)+1);          %Get rid of the edges.
pxm=(px<(info.Nx-(intD/2)));
pr=(pr<(R-1)^2);
pyM=(py>((intD/2)+1));
pym=(py<(info.Ny-(intD/2)));
ii=find((pxM.*pxm.*pr.*pyM.*pym)>0);
px=px(ii);
py=py(ii);



Np=length(ii);
difp=zeros(1,Np);
for np=1:Np
    
    itx=px(np)-intD/2:px(np)+intD/2;
    ity=py(np)-intD/2:py(np)+intD/2;
    simp=sim(ity,itx);
    simop=simo(ity,itx);
    difim=(simp/mean(simp(:))-simop/mean(simop(:))).*circle;
    
    difp(np)=sum(abs(difim(:)));
    
end

maxdifp = max(difp);
participationratio =  sum(difp.^4)/(sum(difp.^2)).^2;

startavalanche = ((maxdifp < cutoffdiff)+(participationratio < cutoffpr)) >0 ;
initial = initial + startavalanche;
%%
nfiles = final-initial+1;

if (nfiles > 0)
    px = zeros(3500,nfiles*350);
    py = zeros(3500,nfiles*350);
    NPF = zeros(1,nfiles*350);
    
    
    
    for ii = initial:final
        fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,ii);
        load(fn,'pxs','pys','Npf');
        MNpf = max(Npf);
        px(1:MNpf,((ii - initial)*350)+1:(ii-initial+1)*350) = pxs(1:MNpf,1:350);
        py(1:MNpf,((ii - initial)*350)+1:(ii-initial+1)*350) = pys(1:MNpf,1:350);
        NPF(((ii - initial)*350)+1:(ii-initial+1)*350) = Npf(1:350);
    end
    
end
%This is how to get initial and final.

% avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En)
% load(avanofile,'avan','navfile');
% changefileindex = find(diff(avan(1,navfile))>1);
% initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
% finalfileindex = navfile(changefileindex);