function [num_part_x num_part_y R_step_first R_step_last] = change_after_rotation(En)
%% Find the change in potential energy after each rotation.

%% Load file that contains both the info about the rotation step and the file numbers
%filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
%avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
%ava_size_file = sprintf('%sAvalanche_size_%i_%i.mat',filedirectory,En,2);
ava_size_file = sprintf('Avalanche_size_%i_%i.mat',En,2);
%load(avanofile);
load(ava_size_file,'Rotation_step','In_imafile','Fn_imafile',...
    'in_trackedfile','Displacement_File_nb','Delta_x','Delta_y','DLength','Initial_CM','Final_CM');
r_step = Rotation_step(1,:);
[R_step_first, ir_first] = unique(r_step,'first'); %First avalanche after rot step n
[R_step_last,  ir_last] = unique(r_step,'last'); %last avalanche after rotation step n
Num_rot = length(R_step_first);  %Number of rotation step with subsequent avalanches
dR_steps = diff(R_step_first);
delta_i = ir_last - ir_first; %number of avalanches after rotation step n.

delta_CM = zeros(2,Num_rot);
delta_r = zeros(2,Num_rot);
    
%% Loop over totation steps
for ir = 1 :Num_rot
    delta_CM(:,ir) = Final_CM(:,ir_last(ir))-Initial_CM(:,ir_first(ir));
    delta_r(:,ir) = [sum(Delta_x(ir_first(ir):ir_last(ir))) sum(Delta_y(ir_first(ir):ir_last(ir)))];
end

%plot(delta_CM(1,:),delta_r(1,:),'.-');    
[d_CMx ,i_unique] = unique(delta_CM(1,:));
num_part_x = fit(d_CMx',delta_r(1,i_unique)','poly1');
[d_CMy ,i_unique] = unique(delta_CM(2,:));
num_part_y = fit(d_CMy',delta_r(2,i_unique)','poly1');
%% This part failed cause it was missing some displacements so the conections wasnt working
% %% Find Avalanches by rotation step and tracked file
% r_step = Rotation_step(1,:);
% [R_step_first, ir_first] = unique(r_step,'first');
% [R_step_last,  ir_last] = unique(r_step,'last');
% Num_rot = length(R_step_first);
% Change_PX = zeros(4000,Num_rot);
% Change_PY = zeros(4000,Num_rot);
% Disk_moved = zeros(1,Num_rot);
% 
% for nr = 1:Num_rot
%     
%     d_file = Displacement_File_nb(_step;
%     Disp_file = sprintf('%sDisplacement_%i.mat',filedirectory,d_file);
%     load(Disp_file,'PX','PY','diskmove');
% 
%     if(ir_first(nr) ~= ir_last(nr))%if info for one rot is in more than one file
%         for i_file = ir_first(nr)+1:ir_last(nr)
%             Disp_file = sprintf('%sDisplacement_%i.mat',filedirectory,i_file);
%             Dpos = load(Disp_file,'PX','PY','diskmove');
%             x1 = PX(:,end);
%             y1 = PY(:,end);
%             x2 = Dpos.PX(:,1);
%             y2 = Dpos.PY(:,1);
%             [adjacentmatrix, trivialbondt1, trivialbondt2] = ...
%                 adjacent(x1,y1,x2,y2,3);
%             PX = PX(trivialbondt1,:);
%             PY = PY(trivialbondt1,:);
%             diskmove = trivialbondt1(diskmove);
%             PX = [PX Dpos.PX(trivialbondt2,:)];
%             PY = [PY Dpos.PY(trivialbondt2,:)];
%             diskmove = [diskmove; trivialbondt2(Dpos.diskmove)];
%         end
%     end
%     
%     Change_PX(diskmove,nr) = PX(diskmove,size(PX,2))-PX(diskmove,1);
%     Change_PY(diskmove,nr) = PY(diskmove,size(PX,2))-PY(diskmove,1);
%     Disk_moved(nr) = length(diskmove);
% end
%     
%     
    
    
    
        

    
