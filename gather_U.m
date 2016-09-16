%% Loads the data from each experiment and saves the change in potential
%energy per avalanche or dissipated energy as well as the avalanche
%duration and the profile change of the pile.



%% Info about file names
%103 104 
%
filenumbers = [15 16 17 18 19 20 21 22 103 104 106]; %Files that contain the info
%filenumbers = [103 104];
Nofiles = length(filenumbers);
kind = 2;


%% Declare variables to save info
CHANGE_U = [];
EXTRA_U = [];
EXTRA_E = [];
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
for ne = 1:Nofiles-1
    En = filenumbers(ne);
    %get info of file to classify the avalanches
    load(sprintf('lost_ava_info_%i.mat',En),'extra_ava_pot','nb_ava_wrong_i','extra_ava_DH');
    
    ava_size_file = sprintf('Avalanche_size_%i_%i.mat',En,kind);
    load(ava_size_file, 'In_imafile','Fn_imafile','Number_Avalanches',...
        'mat_potential','Dheight','Avalanche_duration','diff_Center_mass'...
        ,'Normalized_potential','Num_part_ini','Num_part_end',...
        'in_trackedfile','Nb_boundary','Initial_Angle','Final_Angle',...
        'Avalanche_displacement');
    
    avalanche_nb =  1:Number_Avalanches;
    remain_nb_ava = setdiff(avalanche_nb,nb_ava_wrong_i);
    
    
    EXTRA_U = [EXTRA_U extra_ava_DH(extra_ava_DH<0)];
    EXTRA_E = [EXTRA_E En*ones(size(extra_ava_pot))];
    
    CHANGE_U = [CHANGE_U Dheight(remain_nb_ava)];

    DURATION = [DURATION Avalanche_duration(remain_nb_ava)];
    Pile_Profile_Change = [Pile_Profile_Change ...
        sum(diff_Center_mass(:,remain_nb_ava)>1)/size(diff_Center_mass,1)];
    bottom_of_pile = [bottom_of_pile diff_Center_mass(1:5,remain_nb_ava)];
    U_shape = [U_shape Normalized_potential(:,remain_nb_ava)];
    %Num_part = [Num_part (Num_part_ini+Num_part_end)/2];
    File_num = [File_num in_trackedfile(remain_nb_ava)];
    Experiment_num = [Experiment_num En*ones(size(Dheight(remain_nb_ava)))];
    Num_boundary = [Num_boundary Nb_boundary(:,remain_nb_ava)];
    I_Angle = [I_Angle Initial_Angle(remain_nb_ava)];
    F_Angle = [F_Angle Final_Angle(remain_nb_ava)];
    Angle_change = [Angle_change Final_Angle(remain_nb_ava) - Initial_Angle(remain_nb_ava)];
    DISPLA =[DISPLA Avalanche_displacement(remain_nb_ava)];
    
    
end
% Experiment_num = Experiment_num(2:end);
% i_ch = find(diff(Experiment_num));
% i_ch = [0 i_ch length(Experiment_num)];