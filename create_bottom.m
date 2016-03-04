function nb_files = create_bottom(folder,En)
%% Find the positions of the particles that dont appear in the picture.
%tic;
% Data keeping files
% folder = 1;
% En = 104;
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
% Get files where there were rotation btw them
[r_i_unique, i_rotate,~] = unique(r_i);
i_files_to_rotate = initialfileindex(i_rotate);
f_files_to_rotate = finalfileindex(i_rotate);
%% Saving parameters
N_particles = zeros(1,nb_files);
x_cm = zeros(1,nb_files);
y_cm = zeros(1,nb_files);

 N_rot_particles = zeros(1,nb_files);
 x_from_rot_cm = zeros(1,nb_files);
 y_from_rot_cm = zeros(1,nb_files);
%% Get info about how the center of rotation moves and rotation step
[A, a0_x, a0_y, phi_x] = modify_center_func(folder,En);
phi_y = phi_x + 0.1;
% a0_x = fit_x.a0;
% a0_y = fit_y.a0;
 th_step = 0.0032; %angle per rotation
steps_to_fill = ceil(angle_aperture/th_step);
%% Main loop to complete bottom of the circle of all initialfileindexfiles 
% and get the potential energy before rotation.
for ii = 1:nb_files
    %Extras per image
    x_from_rot = [];
    y_from_rot = [];

    
    file_n = initialfileindex(ii);
    image_fn = sprintf('%s/positions%02d_%05d.mat',filedirectory,En,file_n);
    load(image_fn,'pxs','pys','Npf');
    x_ima = pxs(1:Npf(1),1);
    y_ima = pys(1:Npf(1),1);
    max_before = r_i(ii) - r_i(1); % Max rotation steps before ii
    max_after = r_i(end) - r_i(ii); % Max rotationsteps after ii
    
    d_r = r_i(ii) - r_i_unique;
    ii_new = find(d_r == 0);
    
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
    
    for jj = ii_new-1:-1:i_first_ima
        file_r = f_files_to_rotate(jj);
        image_r = sprintf('%s/positions%02d_%05d.mat',filedirectory,En,file_r);
        load(image_r,'pxs','pys','Npf');
        x = pxs(1:Npf(1),351);
        y = pys(1:Npf(1),351);
        %rotate by last rotation angle to check wich positions to keep
        [~, y_rot] = rotation_composition(x,y,r_i_unique(jj+1),r_i_unique(jj),...
            A,A,phi_x,phi_y,-th_step,a0_x,a0_y);
        x = x(y_rot >= (max(y) - 2*D));
        y = y(y_rot >= (max(y) - 2*D));
        [x_rot, y_rot] = rotation_composition(x,y,r_i(ii),r_i_unique(jj),...
            A,A,phi_x,phi_y,-th_step,a0_x,a0_y);
        x_rot = x_rot(y_rot >= (max(y)-D));
        y_rot = y_rot(y_rot >= (max(y)-D));
        %plot(x_ima,y_ima,'.',x_rot,y_rot,'.');axis('equal');
        %drawnow;
        % Find the ones that are closer than a minimal separation
        adjacentmatrix = adjacent(x_ima, y_ima, x_rot, y_rot, 2/3*D);
        overlap = find(sum(adjacentmatrix));
        i_non_overlap = setdiff((1:length(x_rot)),overlap);
        x_ima = [x_ima ; x_rot(i_non_overlap)];
        y_ima = [y_ima ; y_rot(i_non_overlap)];
        x_from_rot = [x_from_rot ; x_rot(i_non_overlap)];
        y_from_rot = [y_from_rot ; y_rot(i_non_overlap)];
        
    end
    
    
    %Use posterior image
    for jj = ii_new+1:i_last_ima
        file_r = i_files_to_rotate(jj);
        image_r = sprintf('%s/positions%02d_%05d.mat',filedirectory,En,file_r);
        load(image_r,'pxs','pys','Npf');
        x = pxs(1:Npf(1),1);
        y = pys(1:Npf(1),1);
         %rotate by last rotation angle to check wich positions to keep
        [~, y_rot] = rotation_composition(x,y,r_i_unique(jj-1),r_i_unique(jj),...
            A,A,phi_x,phi_y,-th_step,a0_x,a0_y);
        x = x(y_rot >= (max(y) - 2*D));
        y = y(y_rot >= (max(y) - 2*D));
        [x_rot, y_rot] = rotation_composition(x,y,r_i(ii),r_i_unique(jj),...
            A,A,phi_x,phi_y,-th_step,a0_x,a0_y);
        x_rot = x_rot(y_rot >= (max(y)-D));
        y_rot = y_rot(y_rot >= (max(y)-D));
         %plot(x_ima,y_ima,'.',x_rot,y_rot,'.');axis('equal');
        %drawnow;
        % Find the ones that are closer than a minimal separation
        adjacentmatrix = adjacent(x_ima, y_ima, x_rot, y_rot, 2/3*D);
        overlap = find(sum(adjacentmatrix));
        i_non_overlap = setdiff((1:length(x_rot)),overlap);
        x_ima = [x_ima ; x_rot(i_non_overlap)];
        y_ima = [y_ima ; y_rot(i_non_overlap)];
        x_from_rot = [x_from_rot ; x_rot(i_non_overlap)];
        y_from_rot = [y_from_rot ; y_rot(i_non_overlap)];
        
    end
   N_particles(ii) = length(x_ima);
   x_cm(ii) = mean(x_ima);
   y_cm(ii) = mean(y_ima);
   
   N_rot_particles(ii) = length(x_from_rot);
   x_from_rot_cm(ii) = mean(x_from_rot);
   y_from_rot_cm(ii) = mean(y_from_rot);
    
    file_save_i = sprintf('%sComplete_positions_%i.mat',filedirectory,file_n);
    save(file_save_i,'x_ima','y_ima','x_from_rot','y_from_rot');
    
end
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');  
file_save_CM = sprintf('%sCenter_of_Mass_%i.mat',filedirectory, En);
save(file_save_CM,'N_particles','x_cm','y_cm','N_rot_particles',...
    'x_from_rot_cm','y_from_rot_cm','git_version');
            
            
   %toc