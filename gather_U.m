%% Loads the data from each experiment and saves the change in potential 
%energy per avalanche or dissipated energy as well as the avalanche
%duration and the profile change of the pile.


%% Info about file names
% 
filenumbers = [15 16 17 18 19 20 21 22 103 104  106 107]; %Files that contain the info
%filenumbers = [103 104];
Nofiles = length(filenumbers);
kind = 2; 

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
    filename = sprintf('Avalanche_size_%i_%i.mat',En,kind);
    load(filename, 'Dheight','Avalanche_duration','diff_Center_mass'...
        ,'Normalized_potential','Num_part_ini','Num_part_end',...
        'in_trackedfile','Nb_boundary','Initial_Angle','Final_Angle',...
        'Avalanche_displacement');
    CHANGE_U = [CHANGE_U Dheight];
    DURATION = [DURATION Avalanche_duration];
    Pile_Profile_Change = [Pile_Profile_Change sum(diff_Center_mass>1)/size(diff_Center_mass,1)];
    bottom_of_pile = [bottom_of_pile diff_Center_mass(1:5,:)];
    U_shape = [U_shape Normalized_potential];
    Num_part = [Num_part (Num_part_ini+Num_part_end)/2];
    File_num = [File_num in_trackedfile];
    Experiment_num = [Experiment_num En*ones(size(Dheight))];
    Num_boundary = [Num_boundary Nb_boundary];
    I_Angle = [I_Angle Initial_Angle];
    F_Angle = [F_Angle Final_Angle];
    Angle_change = [Angle_change Final_Angle - Initial_Angle];
    DISPLA =[DISPLA Avalanche_displacement];
end
