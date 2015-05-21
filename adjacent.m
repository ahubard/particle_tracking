function [adjacentmatrix, trivialbondt1, trivialbondt2, distancematrix] = adjacent(x1,y1,x2,y2,maxdisp)

%Returns a sparse matrix adjacentmatrix with value ij= 1 if [xi,yi] at time 1 (x1,y1) are with a
%distance <=maxdisp of [xj,yj] at time 2 (x2,y2). distancematrix is the distances between the connected particles.  
%size of x1 and y1 must be (np,1)

%Get size 
Np1 = length(x1);
Np2 = length(x2);

%% Make the bins

initialbinsize = ceil(maxdisp/2)+1;

% Get area where the particles live
xmin = min([x1; x2])-1;
xmax = max([x1; x2])+1;
lx = xmax - xmin;
ymin = min([y1; y2])-1;
ymax = max([y1;y2])+1;
ly = ymax - ymin;

nbinx = max(floor(lx/initialbinsize),1);
nbiny = max(floor(ly/initialbinsize),1);
binsize = max(lx/nbinx , ly/nbiny);


nbinx = nbinx+2;  %Make 1 empty bin at each side of the lx ly area
nbiny = nbiny+2;


%Put particles into the bins. Particle n is in the (ixo(n),iyo(n)) bin
ix1 = ceil((x1-xmin)/binsize)+1;   %+1 to leave the first bin empty
iy1 = ceil((y1-ymin)/binsize)+1;
ix2 = ceil((x2-xmin)/binsize)+1;
iy2 = ceil((y2-ymin)/binsize)+1;
%% Find neighbors of particle np  particlet1 at time one with particlet2 at time t2
%and store its distances

%Create cell array to store which particles are in which bin
boxcontent = cell(nbinx,nbiny);

%Stupid way of getting which particles at time t2 go to which bins. 
for np = 1:Np2
    boxcontent{ix2(np),iy2(np)} = [boxcontent{ix2(np),iy2(np)}, np];
end

particlet1 = [];  
particlet2 = [];
distances = [];

%Find neighbors of particle np.
for np = 1:Np1
    auxiliar = [boxcontent{ix1(np)-1:ix1(np)+1,iy1(np)-1:iy1(np)+1}];
    particlet1 = [particlet1 np*ones(size(auxiliar))];
    particlet2 = [particlet2 auxiliar];
    distances = [distances abs((x1(np)-x2(auxiliar))+1i*(y1(np)-y2(auxiliar)))'];
end

%% Get rid of conections with distances bigger than maxdisp

particlet1 = particlet1(distances<=maxdisp);
particlet2 = particlet2(distances<=maxdisp);
distances = distances(distances<=maxdisp);

%% Create sparse matrices including conectivity and distance. 
distancematrix = sparse(particlet1,particlet2,distances,Np1,Np2);
adjacentmatrix = sparse(particlet1, particlet2, ones(size(particlet1)),Np1,Np2);
trivialbondt1 = find(sum(adjacentmatrix*adjacentmatrix')==1);
%trivialbondt2 = zeros(1, length(trivialbondt1));

[auxtrivial1, auxtrivial2, ~] = find(adjacentmatrix);
[~,itrivial,~] = intersect(auxtrivial1,trivialbondt1);
trivialbondt1 = auxtrivial1(itrivial);
trivialbondt2 = auxtrivial2(itrivial);


% for a=1:length(trivialbondt1) ; 
%     [~, trivialbondt2(a), ~] = find(adjacentmatrix(trivialbondt1(a),:)); 
% end
