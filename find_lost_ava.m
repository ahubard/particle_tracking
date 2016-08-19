%%
function  [nb_ava_wrong_i, new_ava_file1 ,new_ava_file2] = find_lost_ava(En)
% Find out which images have issues ussing the trainning data.
%then check if they cutoff an avalanche or is a new different avalanche.

%% Make classifier
%load trainning set
if(En < 100)
 filenumbers = [15 16 17 18 19 20]; %Files that contain the training info
else
    filenumbers = [103 104 106];
end
Nofiles = length(filenumbers);
kind = 2;
Folder = 1;
filedirectory = sprintf('');

trainning_set = [];
for nex = 1:Nofiles
    %filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',Folder,Folder,filenumbers(ne));
    
    load(sprintf('%slearn_ava_%i.mat',filedirectory, filenumbers(nex)),...
        'ave_ima_diff','part_ratio','part_moved');
    
    trainning_set = [trainning_set [ave_ima_diff;part_ratio;part_moved]];
end

%use trainning set to get a classifier
[trainedClassifier, ~] = ...
    trainClassifier_tofindava(trainning_set);

%% Get data from experiment to fix
%filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
ima_info_file = sprintf('%simages_info_%i.mat',filedirectory,En);
load(ima_info_file,'i_missing_ava','ave_ima_diff','part_ratio',...
    'ima_diff','cutoff1');

%% used classifier to check if there are missing avalanches in such files
i_misshaps = trainedClassifier.predictFcn([ave_ima_diff; part_ratio]);
i_missing_ava = i_missing_ava(i_misshaps > 0);

%% identify to wich image files such misshaps come from
ima_diff_file = sprintf('%sIMA_DIFF_%i.mat',filedirectory,En);
load(ima_diff_file,'file_same_rot1','file_same_rot2','with_undersquare');
file1 = file_same_rot1(i_missing_ava);
file2 = file_same_rot2(i_missing_ava);
N_issues = length(i_missing_ava);
%% find out to which avalanches such files belong to if they do.
%if files with misshaps are in the midle of the avalanches file, get rid of
%such avalanches, if they are in the first or last file then check if it
%cut it out or is a new different avalanche.

%load avalanche image files info
ava_size_file = sprintf('%sAvalanche_size_%i_%i.mat',filedirectory,En,2);
load(ava_size_file,'In_imafile','Fn_imafile','Number_Avalanches',...
    'mat_potential','Avalanche_time');

%initialize variables that will represent

nb_ava_wrong_i = zeros(1,N_issues); %avalanches that something is missing in the middle
new_ava_file1 = zeros(1,N_issues); %file1  avalanche that was missed
new_ava_file2 = zeros(1,N_issues); %file 2 of avalanche that was missed

for ii = 1:N_issues
    i_ini_1 = find(In_imafile <= file1(ii) ,1,'last');
    i_fin_1 = find(Fn_imafile >= file1(ii),1,'first');
    i_ini_2 = find(In_imafile <= file2(ii) ,1,'last');
    i_fin_2 = find(Fn_imafile >= file2(ii),1,'first');
    
    
    if(i_ini_1 == i_fin_1) %file 1 belongs to avalanche i_ini_1
        if (file1(ii) == Fn_imafile(i_fin_1)) %is it at the end?
            num_frames = 350*(Fn_imafile(i_fin_1)-In_imafile(i_fin_1));
            last_ava_frame = find(mat_potential(:,i_fin_1)~=0,1,'last');
            
            if((num_frames - last_ava_frame) < 2) %avalanche i_ni_1 might be wrong
                nb_ava_wrong_i(ii) = i_ini_1;
            else
                if(In_imafile(i_fin_1+1) ~= file2(ii)) %is not the beggining of the next avalanche
                    new_ava_file1(ii) = file1(ii);
                    new_ava_file2(ii) = file2(ii);
                else
                    if(abs(mat_potential(1,i_fin_1+1))>.8) %starts abruptly
                        nb_ava_wrong_i(ii) =  i_fin_1+1;
                    else
                        new_ava_file1(ii) = file1(ii);
                        new_ava_file2(ii) = file2(ii);
                    end
                end
            end
        else %file is in the midlee or beggining and avalanche i_ni_1 is wrong
            nb_ava_wrong_i(ii) = i_ini_1;
        end
    else %file 1 is not part of an avalanche
        if (i_ini_2 == i_fin_2) %if file2 is part of the files of an avalanche
            if(file2(ii) == In_imafile(i_ini_2)) %is it at the beggining?
                if(abs(mat_potential(1,i_ini_2))>0.8)%if it starts abruptly
                    nb_ava_wrong_i(ii) = i_ini_2;%is worng
                else %starts smothly
                    new_ava_file1(ii) = file1(ii); %new avalanche
                    new_ava_file2(ii) = file2(ii);
                end
            else %is in the middle
                nb_ava_wrong_i(ii) = i_ini_2;
            end
        end %file2 is not part of any avalanche files
        new_ava_file1(ii) = file1(ii); %new avalanche
        new_ava_file2(ii) = file2(ii);
    end
    
end

new_ava_file1 = new_ava_file1(new_ava_file1 > 0);
new_ava_file2 = new_ava_file2(new_ava_file2 > 0);
nb_ava_wrong_i = nb_ava_wrong_i(nb_ava_wrong_i > 0);

    
save(sprintf('%slost_ava_info_%i.mat',filedirectory,En),'new_ava_file1',...
    'new_ava_file2','nb_ava_wrong_i')
    
    
    
    
    
    
    
    
    
