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


function [count]=displacement(folder,En,NEn,count,initial,final, alpha, git_version)
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
D = 10;
Ly = 400; %height of image. 
%alpha = 29;  %Angle from the horizontal from experiment.
salpha = sin(alpha*pi/180);
calpha = cos(alpha*pi/180);
%% File to load

filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
%if(En<100)
   % filekernel =sprintf('Tracked_%i%03i',En,NEn);
%else
    filekernel =sprintf('Tracked_%i',NEn);
%end

fnt =sprintf('%s%s.mat',filedirectory,filekernel);

%% 
if(exist(fnt,'file'))
    load(fnt,'PX','PY');
    numframes = size(PX,2); 
    %% Get the particles involved in avalanche.
    iframe1=(PX(:,1)>0);                             %Particles indexes on first frame
    ilastframe=(PX(:,numframes)>0);                  %On last frame
    iframe=find((ilastframe.*iframe1)==1);           %Particles that didnt dissapear in the tracking
    PX = PX(iframe,:);
    PY = Ly - PY(iframe,:); %Mirror the image to get the regular (x,y) axis. 
    %% Get running average position of the particles.
    windowSize = 5;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    x = filter(b,a,PX')';
    y = filter(b,a,PY')';
    
    %% Find the particles that moved.
    totaldelta = abs(x(:,windowSize+1)-x(:,end)+1i*(y(:,windowSize+1)-y(:,end)));
    diskmove = find(totaldelta >0.9);
    %Nummoved = length(diskmove);
    %% Find displacements of such particles btw nframes
    increment=10;
    %dt = increment/2;
    dxraw = (PX(diskmove,1+increment:1:numframes)-PX(diskmove,1:1:numframes-increment))/increment;
    dyraw = (PY(diskmove,1+increment:1:numframes)-PY(diskmove,1:1:numframes-increment))/increment;
    dh = (dyraw*calpha+dxraw*salpha);
    dl = dxraw*calpha-dyraw*salpha;
    drraw2 = dxraw.^2+dyraw.^2;
    dxfil =(x(diskmove,windowSize+increment:1:numframes)-x(diskmove,windowSize:1:numframes-increment))/increment;
    dyfil =(y(diskmove,windowSize+increment:1:numframes)-y(diskmove,windowSize:1:numframes-increment))/increment;
    drfil2 = dxfil.^2+dyfil.^2;
    %dr=sqrt(dr);
    totaltr = ((PX(diskmove,numframes)-PX(diskmove,1)).^2+(PY(diskmove,numframes)-PY(diskmove,1)).^2);
    %totdisplacement=sum(to).*particlesthatmoved)taltr);
    participationratio = sum(totaltr.^4)/(sum(totaltr.^2)^2);
    %% Particles that go over the rim
    Nb_over_boundary = separate_avalanches(filedirectory, En, D, PX,PY);
    
    
    %% Save results
    fnn =sprintf('%sDisplacement_%i.mat',filedirectory,count);
    count = count + 1;
    %fnn =sprintf('Displacement_%i.mat',Count);
    save(fnn,'git_version','windowSize', 'PX','PY','Nb_over_boundary','dh','dl','drraw2','drfil2','diskmove','increment','participationratio','folder','En','NEn','initial','final');
else
    save(sprintf('%sWarning_The_file_%s_does_not_exist.mat',filedirectory,filekernel),'count');
    
end
