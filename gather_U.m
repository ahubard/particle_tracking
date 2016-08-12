%% Loads the data from each experiment and saves the change in potential
%energy per avalanche or dissipated energy as well as the avalanche
%duration and the profile change of the pile.



%% Info about file names
%103 104
%
filenumbers = [15 16 17 18 19 20 21 22]; %Files that contain the info
%filenumbers = [103 104];
Nofiles = length(filenumbers);
kind = 2;

%% Find out which avalanches have issues ussing the trainning data.
%load trainning set
trainning_set = [];
for ne = 1:6
    load(sprintf('learn_ava_%i.mat',filenumbers(ne)),...
        'ave_ima_diff','part_ratio','part_moved');
    
    trainning_set = [trainning_set [ave_ima_diff;part_ratio;part_moved]];
end

%use trainning set to get a classifier
[trainedClassifier, validationAccuracy] = ...
    trainClassifier_tofindava(trainning_set);
%% Declare variables to save info
CHANGE_U = [];
DURATION = [];
Pile_Profile_Change = [];
bottom_of_pile = [];
U_shape = [];
Num_part = [];
File_num =[];
Experiment_num = [];
Num_boundary = [];
I_Angle = [];
F_Angle = [];
Angle_change = [];
DISPLA = [];
%% Main loop over all the files
for ne = 1:Nofiles
    En = filenumbers(ne);
    %get info of file to classify the avalanches
    load(sprintf('images_info_%i.mat',En),'i_missing_ava',...
        'ave_ima_diff','part_ratio');
    %use clasifier to find out if it did missed the avalanche;
    missed_ava = trainedClassifier.predictFcn([ave_ima_diff; part_ratio]);
    %according to the previous classifing reduce the indexes in which and
    %avalanche was missed
    i_missing_ava = i_missing_ava(missed_ava > 0);
    
    %find out to which avalanches this missed files/images correspond
    load(sprintf('IMA_DIFF_%i.mat',En),'file_same_rot1','file_same_rot2');
    %Number of the image file involved in problem.
    f_ima_missing_ava = unique([file_same_rot1(i_missing_ava) file_same_rot2(i_missing_ava)]);
     

    
    ava_size_file = sprintf('Avalanche_size_%i_%i.mat',En,kind);
    load(ava_size_file, 'In_imafile','Fn_imafile','Number_Avalanches',...
        'mat_potential','Dheight','Avalanche_duration','diff_Center_mass'...
        ,'Normalized_potential','Num_part_ini','Num_part_end',...
        'in_trackedfile','Nb_boundary','Initial_Angle','Final_Angle',...
        'Avalanche_displacement');
    
        avalanche_nb =  1:Number_Avalanches;
        %Go over the images of each avalanche and find if the belong to the
        %wrong list. 
    for ifile = 1:Number_Avalanches
         i_ini = find(f_ima_missing_ava >= In_imafile(ifile) ,1,'first');
         i_fin = find(f_ima_missing_ava <= Fn_imafile(ifile),1,'last');
         if(i_ini == i_fin)
            avalanche_nb(ifile) = 0;
         end
    end
    
    %surviving avalanches with no errors from images. 
    
    remain_nb_ava = avalanche_nb(avalanche_nb>0);
    
         
        
        
        
    
    CHANGE_U = [CHANGE_U Dheight(remain_nb_ava)];
    %     DURATION = [DURATION Avalanche_duration]; Pile_Profile_Change =
    %     [Pile_Profile_Change
    %     sum(diff_Center_mass>1)/size(diff_Center_mass,1)]; bottom_of_pile
    %     = [bottom_of_pile diff_Center_mass(1:5,:)]; U_shape = [U_shape
    %     Normalized_potential]; Num_part = [Num_part
    %     (Num_part_ini+Num_part_end)/2]; File_num = [File_num
    %     in_trackedfile];
    Experiment_num = [Experiment_num En*ones(size(Dheight(remain_nb_ava)))];
    %     Num_boundary = [Num_boundary Nb_boundary]; I_Angle = [I_Angle
    %     Initial_Angle]; F_Angle = [F_Angle Final_Angle]; Angle_change =
    %     [Angle_change Final_Angle - Initial_Angle]; DISPLA =[DISPLA
    %     Avalanche_displacement];
    
    
end
Experiment_num = Experiment_num(2:end);
i_ch = find(diff(Experiment_num));
i_ch = [0 i_ch length(Experiment_num)];