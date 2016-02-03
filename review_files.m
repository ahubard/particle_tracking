function files_index = review_files(folder,En)

%% Goes over files and checks if variable IMA exists in it. 

filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
max_num_files = 4001;
if (En == 109)
max_num_files = 3576;
end
variable = sprintf('IMA');

aux_index = zeros(1,max_num_files);

for nf = 1:max_num_files
    if(En > 100)
        filekernel =sprintf('onestep%i_%05i',En,nf);
    else
        filekernel =sprintf('onestep%i%05i',En,nf);
    end
    fno =sprintf('%s%s.mat',filedirectory,filekernel);
    image_in_file = whos(matfile(fno),variable);
    aux_index(nf) = length(image_in_file);
end

files_index = find(aux_index);
fsave = sprintf('%sfileindex_%i.mat',filedirectory,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile,'navfile');

exclude_files = setdiff(files_index,navfile);
save(fsave,'files_index','exclude_files');

%%
% D = 10;
% w = 0.8000;
% Cutoff = 11.5000;
% MinSep = 6.08;




% NF = length(exclude_files);
% Keep = zeros(1,NF);
% Part_rad = zeros(1,NF);
% Max_disp = zeros(1,NF);
% Num_P = zeros(1,NF);
% Std_image = zeros(1,NF);
% for ii = 1 : NF
%     ni = exclude_files(ii);
%     [Keep(ii), ~, ~, ~, ~, ~, Max_disp(ii), Part_rad(ii) ,Num_P(ii), Std_image(ii)] = discriminate(folder,En,ni,D,w,Cutoff,MinSep);
% end
% 
% save(fsave, 'Std_image','Keep','Max_disp','Part_rad','Num_P','-append');




