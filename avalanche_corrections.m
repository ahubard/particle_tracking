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

extra_ava_DH = zeros(1,Nb_new);

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
    
    %% Find particles that didnt move
    [~, trivialbondt1,trivialbondt2] = adjacent(px1,py1,px2,py2,0.9);
    
    %remaining ones
    px1 = px1(setdiff(1:length(px1),trivialbondt1));
    py1 = py1(setdiff(1:length(py1),trivialbondt1));
    px2 = px2(setdiff(1:length(px2),trivialbondt2));
    py2 = py2(setdiff(1:length(py2),trivialbondt2));
    
        
    px1 = px1(py1 < 390);
    py1 = py1(py1 < 390);
    px2 = px2(py2 < 390);
    py2 = py2(py2 < 390);
    
    %% Find the positions in the avalanching part.
    Np1 = length(px1);
    Np2 = length(px2);
    pima_diff_1 = zeros(1,Np1);
    pima_diff_2 = zeros(1,Np2);

    intr = ceil(D/2)+1;
    %
    for np = 1:Np1
        itx1 = px1(np)-intr:px1(np)+intr;
        ity1 = py1(np)-intr:min(py1(np)+intr,400);
        difim1 = imadiff(ity1,itx1);
        pima_diff_1(np) = mean(abs(difim1(:)));
    end
    
    
    for np = 1:Np2
        itx2 = px2(np)-intr:px2(np)+intr;
        ity2 = py2(np)-intr:min(py2(np)+intr,400);
        difim2 = imadiff(ity2,itx2);
        pima_diff_2(np) = mean(abs(difim2(:)));
    end
    

    
    cutoffima = 0.075;
    moved1 = find(pima_diff_1 > cutoffima);
    moved2 = find(pima_diff_2 > cutoffima);
   
    %check if they really moved
    if(~isempty(moved1) && ~isempty(moved2))
        
        %find the one with more moved particles above the cutoff and keep that
        %amount of particles.
        if(length(moved1) > length(moved2))
            [~, i_sort] = sort(pima_diff_2,'descend');
            moved2 = i_sort(1:length(moved1));
        else
            [~, i_sort] = sort(pima_diff_1,'descend');
            moved1 = i_sort(1:length(moved2));
        end
        
        px1 = px1(moved1);
        py1 = py1(moved1);
        px2 = px2(moved2);
        py2 = py2(moved2);
        
        %change refernce system to gravitational pg is position in the gravity
        %direction and pl in the perpendicular to gravity.
        
        
        [pl1, pg1] = rot_me(alpha,px1,py1);
        [pl2, pg2] = rot_me(alpha,px2,py2);
        
        
        nb1 = length(pl1);
        nb2 = length(pl2);
        
        if(~isempty(px1) && ~isempty(px2))
            [~,tb1,tb2,~] = adjacent(pl1,pg1,pl2,pg2,D/2);
            pl1_tb = pl1(tb1);
            pg1_tb = pg1(tb1);
            pl2_tb = pl2(tb2);
            pg2_tb = pg2(tb2);
            
            if ((length(tb1) < nb1) && (length(tb2) < nb2))
                pl1_remain = pl1(setdiff((1:nb1),tb1));
                pg1_remain = pg1(setdiff((1:nb1),tb1));
                pl2_remain = pl2(setdiff((1:nb1),tb1));
                pg2_remain = pg2(setdiff((1:nb1),tb1));

                [~,~,~,distmat] = ...
                    adjacent(pl1_remain,pg1_remain,pl2_remain,pg2_remain,3*D);

                distmat(distmat==0) = inf;
                [assignment, ~] = assignmentoptimal(distmat);

                i_stay = find(assignment>0);
                pl1_aux= pl1_remain(i_stay);
                pl2_aux = pl2_remain(assignment(i_stay));
                pg1_aux = pg1_remain(i_stay);
                pg2_aux = pg2_remain(assignment(i_stay));
            else
                pg1_aux =[];
                pg2_aux =[];

            end
                
        end
        
        pg2 = [pg2_tb(:); pg2_aux(:)];
        pg1 = [pg1_tb(:); pg1_aux(:)];

        
        
        DH = pg2 - pg1;
        extra_ava_DH(in) = sum(DH);
    else
        extra_ava_DH(in) = 0;
    end
    
    
%     
%     Np = length(moved);
%     
%     p_moved = zeros(1,Np);
%     intr = ceil(D/2);
%     
%     for np=1:Np
%         itx1 = PX(1,np)-intr:PX(1,np)+intr;
%         ity1 = PY(1,np)-intr:min(PY(1,np)+intr,400);
%         itx2 = PX(2,np)-intr:PX(2,np)+intr;
%         ity2 = PY(2,np)-intr:min(PY(2,np)+intr,400);
%         difim1 = imadiff(ity1,itx1);
%         difim2 = imadiff(ity2,itx2);
%         p_moved(np) = mean(abs(difim1(:)))+mean(abs(difim2(:)));
%     end
%     
%     %check which ones moved from the image
%     cutoffima = 0.1;
%     moved = find(p_moved > cutoffima);
%     PX = PX(:,moved);
%     PY = PY(:,moved);
%     DX = PX(2,:) - PX(1,:);
%     DY = PY(2,:) - PY(1,:);
%     
%     
%    
end

save(sprintf('%s%s',filedirectory, lost_ava_info),...
    'extra_ava_DH','new_ava_file1','new_ava_file2','nb_ava_wrong_i','-append');
