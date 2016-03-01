%% Find the positions of the particles that dont appear in the picture.
tic;
% Data keeping files
folder = 1;
En = 104;
D = 10;
angle_aperture = 2.35; %Around how much is missing of angle.
%% Get info about the file numbers and rotation steps
filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile);
% navfile has the numbers of files that contain avalanches. avan(1,:)
% indicates if the files are continuos in time or broke. If
% avan(1,n)=avan(1,n+1) then files n and n+1 belong to a continuos set of time.
changefileindex = find(diff(avan(1,navfile))>1);  
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);%first file of set
finalfileindex = navfile(changefileindex);%last file of set. 
% Task is to recreate the bottom of the first image of files
% initialfileindex.
nb_files = length(initialfileindex);
% Between finalfileindex(n) and initialfileindex(n+1) a rotation ocurred.
% The rotation steps t of the files intialfileindex and finalfile index are
% given by: 
r_i = avan(2,initialfileindex);
r_f = avan(2,finalfileindex);
delta_r = r_i(2:end) - r_f;  %number of rotation steps between avalanches
%% Get info about how the center of rotation moves and rotation step size
center_file = sprintf('%sCenter_%i.mat',filedirectory,En);
load(center_file);
a0_x = fit_x.a0;
a0_y = fit_y.a0;
th_step = fit_x.w*.7+fit_y.w*.3; %angle per rotation
steps_to_fill = angle_aperture/th_step;
%% Main loop to complete bottom of the circle of all initialfileindexfiles 
% and get the potential energy before rotation.
for ii = 1:nb_files
    file_n = initialfileindex(ii);
    image_fn = sprintf('%s/positions%02d_%05d.mat',filedirectory,file_n);
    load(image_fn,'pxs','pys','Npf');
    x_ima = pxs(1:Npf(1),1);
    y_ima = pys(1:Npf(1),1);
    max_before = r_i(ii) - r_i(1); %Max rotation steps before ii
    max_after = r_i(end) - r_i(ii); %Max rotationsteps after ii
    d_r = r_i(ii) - r_i(1:end);
    %Find how many rotation steps before and after must be used depending on the
    %current image rotation step.
    if(max_before <= floor(steps_to_fill/2))
        t_r_before = max_before;
    elseif (max_after <= floor(steps_to_fill/2))
        t_r_before = steps_to_fill - max_after;
    else
        t_r_before = ceil(steps_to_fill/2);
    end
    t_r_after = steps_to_fill - t_r_before;
    
    % Get the first and last image to be used
    i_first_ima = find(d_r >= t_r_before,1,'last');
    i_last_ima = find(-d_r >= t_r_after,1,'first');
    
    
    % Use previous images
    if(i_first_ima)
        for jj = ii-1:-1:i_first_ima
            file_r = finalfileindex(jj);
            image_r = sprintf('%s/positions%02d_%05d.mat',filedirectory,file_r);
            load(image_fn,'pxs','pys','Npf');
            x = pxs(1:Npf(1),1);
            y = pys(1:Npf(1),1);
            %rotate by last rotation angle to check wich positions to keep
            [~, y_rot] = rotation_composition(x,y,r_i(jj+1),r_i(jj),...
                A_x,A_y,phi_x,-th_step,a0_x,a0_y);
            x = x(y_rot >= (max(y) - D));
            y = y(y_rot >= (max(y) - D));
            [x_rot, y_rot] = rotation_composition(x,y,r_i(ii),r_i(jj),...
                A_x,A_y,phi_x,-th_step,a0_x,a0_y);
            % Find the ones that overlap or come from both files
            [~, trivialbondt1,trivialbondt2] = adjacent(x_ima,y_ima,x_rot,y_rot,D/2);
            x_both = (x_rot(trivialbondt2) + x_ima(trivialbondt1))/2;
            y_both = (y_rot(trivialbondt2) + y_ima(trivialbondt1))/2;
            %Find the particles that are not in both sets. 
            [x_ima_e, ix] = setdiff(x_ima,x_ima(trivialbondt1));
            y_ima_e = y_ima(ix);
            [x_rot_e, ix] = setdiff(x_rot,x_rot(trivialbondt2));
            y_rot_e = rot_y(ix);
            %Positions from all previous rotations
            x_ima = [x_both; x_ima_e; x_rot_e];
            y_ima = [y_both; y_ima_e; y_rot_e];
            %x_ima(length(x_ima)+1:length(x_ima)+length(x)) = x_rot;
            %y_ima(length(y_ima)+1:length(y_ima)+length(y)) = y_rot;
        end
    end
    
    %Use posterior images
    if(i_last_ima)
        for jj = ii+1:i_last_ima
            file_r = initialfileindex(jj);
            image_r = sprintf('%s/positions%02d_%05d.mat',filedirectory,file_r);
            load(image_fn,'pxs','pys','Npf');
            x = pxs(1:Npf(1),1);
            y = pys(1:Npf(1),1);
            %rotate by last rotation angle to check wich positions to keep
            [~, y_rot] = rotation_composition(x,y,r_i(jj-1),r_i(jj),...
                A_x,A_y,phi_x,-th_step,a0_x,a0_y);
            x = x(y_rot >= (max(y)-10));
            y = y(y_rot >= (max(y)-10));
            [x_rot, y_rot] = rotation_composition(x,y,r_i(ii),r_i(jj),...
                A_x,A_y,phi_x,-th_step,a0_x,a0_y);
            % Find the ones that overlap or come from both files
            [~, trivialbondt1,trivialbondt2] = adjacent(x_ima,y_ima,x_rot,y_rot,D/2);
            x_both = (x_rot(trivialbondt2) + x_ima(trivialbondt1))/2;
            y_both = (y_rot(trivialbondt2) + y_ima(trivialbondt1))/2;
            %Find the particles that are not in both sets. 
            [x_ima_e, ix] = setdiff(x_ima,x_ima(trivialbondt1));
            y_ima_e = y_ima(ix);
            [x_rot_e, ix] = setdiff(x_rot,x_rot(trivialbondt2));
            y_rot_e = rot_y(ix);
            %Positions from all previous rotations
            x_ima = [x_both; x_ima_e; x_rot_e];
            y_ima = [y_both; y_ima_e; y_rot_e];
            %x_ima(length(x_ima)+1:length(x_ima)+length(x)) = x_rot;
            %y_ima(length(y_ima)+1:length(y_ima)+length(y)) = y_rot;
        end
    end
    
    file_save = sprintf('%sComplete_positions_%i.mat',filedirectory,file_n);
    save(file_save,'x_ima','y_ima');
    
end
    
            
            
   toc