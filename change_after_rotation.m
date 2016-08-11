function [Num_rot] = change_after_rotation(folder,En)
%% Find the change in potential energy after each rotation.

%% Load file that contains both the info about the rotation step and the file numbers
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
ava_size_file = sprintf('%sAvalanche_size_%i_%i.mat',filedirectory,En,2);
file_for_center = sprintf('%sCenter_%i.mat',filedirectory, En);
%ava_size_file = sprintf('Avalanche_size_%i_%i.mat',En,2);
%file_for_center = sprintf('Center_%i.mat',En);
load(file_for_center,'A_x','A_y','phi_x','phi_y','frequency','a0_x','a0_y');
load(avanofile,'alpha');
alpha = alpha*pi/180;
load(ava_size_file,'Rotation_step','In_imafile','Fn_imafile',...
    'in_trackedfile','Displacement_File_nb','Delta_x','Delta_y','DLength',...
'Initial_CM','Final_CM','Dheight');
r_step = Rotation_step(1,:);
[R_step_first, ir_first] = unique(r_step,'first'); %First avalanche after rot step n
[R_step_last,  ir_last] = unique(r_step,'last'); %last avalanche after rotation step n
Num_rot = length(R_step_first);  %Number of rotation step with subsequent avalanches
dR_steps = diff(R_step_first);
delta_i = ir_last - ir_first; %number of avalanches after rotation step n.

delta_CM = zeros(2,Num_rot);
delta_r = zeros(2,Num_rot);
delta_g = zeros(2,Num_rot);
x_cm = zeros(1,Num_rot);
y_cm = zeros(1,Num_rot);

%% Loop over rotation steps
for ir = 1 :Num_rot
    delta_CM(:,ir) = Final_CM(:,ir_last(ir))-Initial_CM(:,ir_first(ir));
    delta_r(:,ir) = [sum(Delta_x(ir_first(ir):ir_last(ir))) sum(Delta_y(ir_first(ir):ir_last(ir)))];
    delta_g(:,ir) = [sum(DLength(ir_first(ir):ir_last(ir))) sum(Dheight(ir_first(ir):ir_last(ir)))];
    if (ir < Num_rot)
    [x_cm(ir), y_cm(ir)] = ...
        rotation_composition(Final_CM(1,ir_last(ir)),Final_CM(2,ir_last(ir))...
        ,R_step_first(ir),R_step_first(ir+1),A_x,A_y,phi_x,phi_y, -frequency,a0_x,a0_y);
    end
end

%plot(delta_CM(1,:),delta_r(1,:),'.-');    
[d_CMx ,i_unique] = unique(delta_CM(1,:));
num_part_x = fit(d_CMx',delta_r(1,i_unique)','poly1');
[d_CMy ,i_unique] = unique(delta_CM(2,:));
num_part_y = fit(d_CMy',delta_r(2,i_unique)','poly1');

%% Find Avalanches by rotation step and tracked file
% This part failed cause it was missing some displacements so the conections wasnt working
r_step = Rotation_step(1,:);
[R_step_first, ir_first] = unique(r_step,'first');
[R_step_last,  ir_last] = unique(r_step,'last');
Num_rot = length(R_step_first);
Change_PX = zeros(4000,Num_rot);
Change_PY = zeros(4000,Num_rot);
Disk_moved = zeros(1,Num_rot);
 
%save the displacement of the particles that weren't conected but moved.
DXS = zeros(1000,Num_rot);
DYS = zeros(1000,Num_rot);
fail_to_conect = zeros(1,Num_rot);

for nr = 1:Num_rot
    l_notfound = 0;
    b_file = Displacement_File_nb(ir_first(nr));
    Disp_file = sprintf('%sDisplacement_%i.mat',filedirectory,b_file);
    load(Disp_file,'PX','PY','diskmove')

    if(ir_first(nr) ~= ir_last(nr))%if info for one rot is in more than one file
        for i_file = ir_first(nr)+1:ir_last(nr)
            d_file = Displacement_File_nb(i_file);
            if(d_file ~= b_file)
                Disp_file = sprintf('%sDisplacement_%i.mat',filedirectory,d_file);
                Dpos = load(Disp_file,'PX','PY','diskmove');
                x1 = PX(:,end);
                y1 = PY(:,end);
                x2 = Dpos.PX(:,1);
                y2 = Dpos.PY(:,1);
                [adjacentmatrix, trivialbondt1, trivialbondt2] = ...
                    adjacent(x1,y1,x2,y2,4);
                
                
                [~,tb1_diskmove] = intersect(trivialbondt1,diskmove);
                notfound1 = setdiff(diskmove,trivialbondt1);
                DXS(l_notfound+(1:length(notfound1)),nr) = ...
                    PX(notfound1,end) - PX(notfound1,1);
                DYS(l_notfound+(1:length(notfound1)),nr) = ...
                    PY(notfound1,end) - PY(notfound1,1);
                l_notfound = l_notfound + length(notfound1);
                
                PX = PX(trivialbondt1,:);
                PY = PY(trivialbondt1,:);
                
                
                [~,tb2_diskmove] = intersect(trivialbondt2,Dpos.diskmove);
                diskmove = union(tb1_diskmove,tb2_diskmove);
                notfound2 = setdiff(Dpos.diskmove,trivialbondt2);
                DXS(l_notfound+(1:length(notfound2)),nr) = ...
                    Dpos.PX(notfound2,end) - Dpos.PX(notfound2,1);
                DYS(l_notfound+(1:length(notfound2)),nr) = ...
                    Dpos.PY(notfound2,end) - Dpos.PY(notfound2,1);
                l_notfound = l_notfound + length(notfound2);
                
                PX = [PX Dpos.PX(trivialbondt2,:)];
                PY = [PY Dpos.PY(trivialbondt2,:)];
                
            end
            b_file = d_file;
        end
    end
    
    Change_PX(diskmove,nr) = PX(diskmove,size(PX,2))-PX(diskmove,1);
    Change_PY(diskmove,nr) = PY(diskmove,size(PX,2))-PY(diskmove,1);
    Disk_moved(nr) = length(diskmove);
    fail_to_conect(nr) = l_notfound;
end
    
%in the camera frame reference.
    X2_displacement = sum(Change_PX.^2) + sum(DXS.^2);
    Y2_displacement = sum(Change_PY.^2) + sum(DYS.^2);
    
% in the gravitational frame reference;
[Change_orthog, Change_g] = rot_me(alpha,Change_PX(:),Change_PY(:));
[Dorthog, Dg] = rot_me(alpha,DXS(:),DYS(:));

Orthog_displacement = sum(reshape(Change_orthog,size(Change_PX))) + ...
    sum(reshape(Dorthog, size(DXS)));
G_displacement = sum(reshape(Change_g,size(Change_PY))) + ...
    sum(reshape(Dg,size(DYS)));
    
save_file = sprintf('%sAva_after_rot_%i.mat',filedirectory,En);
    

save(save_file,'X2_displacement','Y2_displacement','Orthog_displacement',...
'G_displacement','R_step_first','ir_first','ir_last','delta_i','delta_CM',...
'delta_r','x_cm','y_cm','fail_to_conect','delta_g');
        

    
