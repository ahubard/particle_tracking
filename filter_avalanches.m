function remain_nb_ava = filter_avalanches(folder,En)
%% Go over the files with issues and keep only the avalanches that are far
% from problems, such as the difference with the file before in the same
% rotation is big so there was some movement moved or it doesnt start and
% ends close to zero but there is an abrupt change.

%% Info of experiment

cutoff1 = 7e5; %cutoff for ima_diff
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);

%% Get the image difference of files within the same rotation
ava_size_file = sprintf('%sAvalanche_size_%i_%i.mat',filedirectory,En,2);
load(ava_size_file,'In_imafile','Fn_imafile','Number_Avalanches',...
    'mat_potential');
ima_diff_file = sprintf('%sIMA_DIFF_%i.mat',filedirectory,En);
if (exist(ima_diff_file))
    load(ima_diff_file,'ima_diff','file_same_rot1','file_same_rot2','with_undersquare');
    i_missing_ava = find(ima_diff > cutoff1);
    nb_miss = length(i_missing_ava);
    
    %Look into the images and apply two extra criterions:
    %1) the diference of the images normalized by the average value of that
    %image
    %2) Check if a simple track conects the positions of the particles
    
    ave_ima_diff = zeros(1,nb_miss);
    part_ratio = zeros(1,nb_miss);
    if(sum(with_undersquare) == 0)
        filebase1 = sprintf('%sonestep%i',filedirectory,En);
    else
        filebase1 = sprintf('%sonestep%i_',filedirectory,En);
    end
    filebase2 = sprintf('%spositions%i',filedirectory,En);
    
 for nf = 1:nb_miss
     %get and compare images
    fn1 = file_same_rot1(i_missing_ava(nf));
    fn2 = file_same_rot2(i_missing_ava(nf));
    I1 = load(sprintf('%s%05i.mat',filebase1,fn1),'IMA');
    I2 = load(sprintf('%s%05i.mat',filebase1,fn2),'IMA');
    im1 = I1.IMA(:,:,351);
    im2 = I2.IMA(:,:,1);
    ave_ima_diff(nf) = sum((im1(:)/mean(im1(:))-im2(:)/mean(im2(:))).^2);
    part_ratio(nf) = sum((im1(:)/mean(im1(:))-im2(:)/mean(im2(:))).^4)/...
        sum((im1(:)/mean(im1(:))-im2(:)/mean(im2(:))).^2)^2;
    
    
%     %get and compare positions
%     P1 = load(sprintf('%s%05i.mat',filebase2,fn1),'pxs','pys','Npf');
%     P2 = load(sprintf('%s%05i.mat',filebase2,fn2),'pxs','pys','Npf');
%     x1 = P1.pxs(1:P1.Npf(351),351);
%     y1 = P1.pys(1:P1.Npf(351),351);
%     x2 = P2.pxs(1:P2.Npf(1),1);
%     y2 = P2.pys(1:P2.Npf(1),1);
%     [adjacentmatrix, trivialbondt1, trivialbondt2, distancematrix] = ...
%         adjacent(x1,y1,x2,y2,4)

 
end
    
    %Number of the image file involved in problem.
    f_ima_missing_ava = unique([file_same_rot1(i_missing_ava) file_same_rot2(i_missing_ava)]);
   
    
    %% Go over all the displacement files to assign a displacement file to each
    %ima file.
    
    
    avalanche_nb = zeros(1,max(f_ima_missing_ava));
    
    
    
    %Get rid of the avalanches avalanche_nb(f_ima_missing_Ava) cause there was
    %an error there.
    
    remain_nb_ava = setdiff(1:Number_Avalanches, avalanche_nb(f_ima_missing_ava));
else
    remain_nb_ava = 1:Number_Avalanches;
end
%% Look at the start and end points of the remaining avalanches to make sure
% they werent cut abruptly or started abruptly
%Keep the avalanches that didnt have problems in the previous module
AVA = mat_potential(:,remain_nb_ava);

%Find the first and last non zero value of each avalanche.
Nb_av = size(AVA,2);
first_val = zeros(1,Nb_av);
last_val = zeros(1,Nb_av);
for na = 1:Nb_av;
    i_first = find(AVA(:,na) ~= 0,1,'first');
    first_val(na) = AVA(i_first,na);
    i_last = find(AVA(:,na) ~= 0,1,'last');
    last_val(na) = AVA(i_last,na);
end

%% Keep only the ones that start and end smoothly
cutoff2 = 0.9;
remain_nb_ava = remain_nb_ava((abs(first_val) < cutoff2)  & (abs(last_val) < cutoff2));

save(ava_size_file,'remain_nb_ava','-append')




