function  Nb_new = avalanche_corrections(folder,En)
%% Files to load
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
lost_ava_info = sprintf('lost_ava_info_%i.mat',En);
bgfile = sprintf('%sback%i.mat',filedirectory,En);

if (En > 100)
    D = 10;
    w = 0.8;
else
    D = 9;
    w = 0.7;
end
%% load find_lost ava, to see whihc avalanches are wrong and which ones must be measured.
load (avanofile,'alpha');
alpha = alpha*pi/180;

load(lost_ava_info);
Nb_new = length(new_ava_file1);

extra_ava_pot = zeros(1,Nb_new);

for in = 1:Nb_new
    clear('PX','PY');
    % Get the imafile
    if (En < 100)
        ima_file_1 = sprintf('%sonestep%i%05i.mat',filedirectory,En,new_ava_file1(in));
    else
         ima_file_1 = sprintf('%sonestep%i_%05i.mat',filedirectory,En,new_ava_file1(in));
    end
    load(ima_file_1,'IMA');
    im1 = IMA(:,:,351);
    %find the missing particles
    [px1, py1, Nbp, sim1] = particles_in_frame(im1,filedirectory,En,D,w);
    
    if(En < 100)
        ima_file_2 = sprintf('%sonestep%i%05i.mat',filedirectory,En,new_ava_file2(in));
    else
       ima_file_2 = sprintf('%sonestep%i_%05i.mat',filedirectory,En,new_ava_file2(in));
    end
        load(ima_file_2,'IMA');
    im2 = IMA(:,:,1);
    %find the missing particles
    [px2, py2, Nbp2, sim2] = particles_in_frame(im2,filedirectory,En,D,w);
    
    imadiff = sim2./mean(sim2(:))-sim1./mean(sim1(:));
    
    
    [~, trivialbondt1,trivialbondt2] = adjacent(px1,py1,px2,py2,D/2);
    
    npf = length(trivialbondt1);
    %trivial ones
    PX(1,:) = px1(trivialbondt1);
    PX(2,:) = px2(trivialbondt2);
    PY(1,:) = py1(trivialbondt1);
    PY(2,:) = py2(trivialbondt2);
    
    %remaining ones
    px1 = px1(setdiff(1:length(px1),trivialbondt1));
    py1 = py1(setdiff(1:length(py1),trivialbondt1));
    px2 = px2(setdiff(1:length(px2),trivialbondt2));
    py2 = py2(setdiff(1:length(py2),trivialbondt2));
    
    px1 = px1(py1 < 370);
    py1 = py1(py1 < 370);
    px2 = px2(py2 < 370);
    py2 = py2(py2 < 370);
    
    
    if(~isempty(px1) && ~isempty(px2))
        [~, trivialbondt1,trivialbondt2,distmat] = adjacent(px1,py1,px2,py2,2.5*D);
        
        distmat(distmat==0) = inf;
        [assignment, ~] = assignmentoptimal(distmat);
        
        i_stay = find(assignment>0);
        np_aux = length(i_stay);
        PX(1,npf+(1:np_aux)) = px1(i_stay);
        PX(2,npf+(1:np_aux)) = px2(assignment(i_stay));
        PY(1,npf+(1:np_aux)) = py1(i_stay);
        PY(2,npf+(1:np_aux)) = py2(assignment(i_stay));
    end
    
    %keep only the ones that moved.
    moved = find((abs(PX(2,:)-PX(1,:))+abs(PY(2,:)-PY(1,:))) > 0);
    PX = PX(:,moved);
    PY = PY(:,moved);
    
    Np = length(moved);
    
    p_moved = zeros(1,Np);
    intr = ceil(D/2);
    
    for np=1:Np
        itx1 = PX(1,np)-intr:PX(1,np)+intr;
        ity1 = PY(1,np)-intr:min(PY(1,np)+intr,400);
        itx2 = PX(2,np)-intr:PX(2,np)+intr;
        ity2 = PY(2,np)-intr:min(PY(2,np)+intr,400);
        difim1 = imadiff(ity1,itx1);
        difim2 = imadiff(ity2,itx2);
        p_moved(np) = mean(abs(difim1(:)))+mean(abs(difim2(:)));
    end
    
    %check which ones moved from the image
    cutoffima = 0.1;
    moved = find(p_moved > cutoffima);
    PX = PX(:,moved);
    PY = PY(:,moved);
    DX = PX(2,:) - PX(1,:);
    DY = PY(2,:) - PY(1,:);
    
    
    [~, change_g] = rot_me(alpha,DX,DY);
    extra_ava_pot(in) = sum(change_g);
end

save(sprintf('%s%s',filedirectory, lost_ava_info),'extra_ava_pot',...
    'new_ava_file1','new_ava_file2','nb_ava_wrong_i');
