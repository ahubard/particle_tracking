function nb_files = create_bottom(folder,En)
%% Find the positions of the particles that dont appear in the picture.
%tic;
% Data keeping files
% folder = 1;
% En = 104;
D = 10;
angle_aperture = 2.5; %Around how much is missing of angle.
m = 0.2;
move_cutoff = 4;
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
    Track_file = sprintf('%sDisplacement_%i.mat',filedirectory, ii);
    load(Track_file,'initial','final','diskmove','PX','PY');
    
    if(file_n ~= initial)
        save(sprintf('%sError_%i.mat',filedirectory,ii),'file_n','initial');
        error('Somethings not right');
    end
    
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
    mtheta = -tan(t_r_before*th_step);
    x_rot = [];
    y_rot = [];
    
    % Use previous images
    tic
    for jj = ii_new-1:-1:i_first_ima
        file_r = f_files_to_rotate(jj);
        Track_file = sprintf('%sDisplacement_%i.mat',filedirectory, i_rotate(jj));
        load(Track_file,'initial','final','diskmove','PX','PY');
        
        if(file_r ~= final)
            save(sprintf('%sError_%i.mat',filedirectory,ii),'file_n','initial','final','jj');
            error('Somethings not right');
        end
        PY = 400 - PY;
        position_change = sqrt((PX(:,end)-PX(:,1)).^2+(PY(:,end)-PY(:,1)).^2);
        particles_move = find(position_change >= move_cutoff);
        
        
        image_r = sprintf('%s/positions%02d_%05d.mat',filedirectory,En,file_r);
        load(image_r,'pxs','pys','Npf');
        [~, t] = max(Npf);
        x = pxs(1:Npf(t),t);
        y = pys(1:Npf(t),t);
        r = sqrt((x-a0_x).^2+(y-a0_y).^2);
        
        %Get the angle between two rotations and make sure to include at
        %least the particles between. 
        alpha_btw_rot = (r_i_unique(jj+1) - r_i_unique(jj))*th_step;
        needed = find((x < a0_x)  &  (y > 400 - r*sin(alpha_btw_rot -D)));
        
        % Keep particles that didnt move much in previous frames
        adjacentmatrix = adjacent(x, y, PX(particles_move,end), PY(particles_move,end), D/2+1);
        particles_move = find(sum(adjacentmatrix,2));
        p_to_rotate = setdiff(1:length(x),particles_move);
        p_to_rotate = union(p_to_rotate,needed);
        
        x_aux = x(p_to_rotate);
        y_aux = y(p_to_rotate);
        r_aux = r(p_to_rotate);
        
        y_cutoff = max(y)-D+m*(x_aux-a0_x);
        y_cutoff(y_cutoff < max(y) - 10*D) = max(y)-12*D;
        y_cutoff = min(y_cutoff, 400-r_aux*sin(alpha_btw_rot)-D);
        
        
        x = x_aux(x_aux < a0_x-D  & y_aux > y_cutoff);
        y = y_aux(x_aux < a0_x-D  & y_aux > y_cutoff);
        
        %rotate positions to meet ima positions
        [x_rot, y_rot] = rotation_composition(x,y,r_i(ii),r_i_unique(jj),...
            A,A,phi_x,phi_y,th_step,a0_x,a0_y);
        x_rot = x_rot(y_rot >= (max(y)-3*D));
        y_rot = y_rot(y_rot >= (max(y)-3*D));
%          plot(x_ima,y_ima,'.',x_rot,y_rot,'.');axis('equal');
%          drawnow;
        [x_ima, y_ima, x_non_ima, y_non_ima, o_s] = merge_positions(x_ima,y_ima,x_rot,y_rot);
        if (o_s(1) == 0 || o_s(2) == 0)
            save(sprintf('%sCheck_bottom_%i.mat',filedirectory, file_n),'file_r')
        end
        x_from_rot = [x_from_rot(:); x_non_ima(:)];
        y_from_rot = [y_from_rot(:); y_non_ima(:)];
        
    end
    [xo, ixo] = max([a0_x; x_rot(y_rot < 400)]);
    
    yo = [a0_y; y_rot(y_rot < 400)];
    yo = yo(ixo);
    %Use posterior image
    for jj = ii_new+1:i_last_ima
        file_r = i_files_to_rotate(jj);
        Track_file = sprintf('%sDisplacement_%i.mat',filedirectory, i_rotate(jj));
        load(Track_file,'initial','final','diskmove','PX','PY');
        
        if(file_r ~= initial)
            save(sprintf('%sError_%i.mat',filedirectory,ii),'file_n','initial','final','jj');
            error('Somethings not right');
        end
        PY = 400 - PY;
        position_change = sqrt((PX(:,end)-PX(:,1)).^2+(PY(:,end)-PY(:,1)).^2);
        particles_move = find(position_change >= move_cutoff);
        
        image_r = sprintf('%s/positions%02d_%05d.mat',filedirectory,En,file_r);
        load(image_r,'pxs','pys','Npf');
        [~, t] = max(Npf);
        x = pxs(1:Npf(t),t);
        y = pys(1:Npf(t),t);
        r = sqrt((x-a0_x).^2+(y-a0_y).^2);
        
        %Get the angle between two rotations and make sure to include at
        %least the particles between. 
        alpha_btw_rot = (r_i_unique(jj) - r_i_unique(jj-1))*th_step;
        needed = find((x > a0_x)  &  (y > 400 - r*sin(alpha_btw_rot - D)));
        
        
        % Keep particles that didnt move much in previous frames
        adjacentmatrix = adjacent(x, y, PX(particles_move,1), PY(particles_move,1), D/2 + 1);
        particles_move = find(sum(adjacentmatrix,2));
        p_to_rotate = setdiff(1:length(x),particles_move);
        p_to_rotate = union(p_to_rotate, needed);
        
        x_aux = x(p_to_rotate);
        y_aux = y(p_to_rotate);
        r_aux = r(p_to_rotate);
        
        y_cutoff = max(y)-D-m*(x_aux-a0_x);
        %y_cutoff(y_cutoff < max(y) - 10*D) = max(y)-12*D;
        y_cutoff = min(y_cutoff, 400-r_aux*sin(alpha_btw_rot)-D);
        
        
        x = x_aux(x_aux > (a0_x+D)  & y_aux > y_cutoff);
        y = y_aux(x_aux > (a0_x+D)  & y_aux > y_cutoff);
        
        % x = x_aux(x_aux < a0_x-D  & y_aux > max(y)-1.5*D+m*(x_aux-a0_x) & y_aux > max(y)-15*D);
        
        %rotate positions to meet ima positions
        [x_aux, y_aux] = rotation_composition(x,y,r_i(ii),r_i_unique(jj),...
            A,A,phi_x,phi_y,th_step,a0_x,a0_y);
        x_rot = x_aux(x_aux >= 1/mtheta*(y_aux-yo)+xo-5*D  & y_aux >(max(y)-3*D));
        y_rot = y_aux(x_aux >= 1/mtheta*(y_aux-yo)+xo-5*D  & y_aux >(max(y)-3*D));
%          plot(x_ima,y_ima,'.',x_rot,y_rot,'.');axis('equal');
%          drawnow;
        [x_ima, y_ima, x_non_ima, y_non_ima, o_s] = merge_positions(x_ima,y_ima,x_rot,y_rot);
        if (o_s(1) == 0 || o_s(2) == 0)
            save(sprintf('%sCheck_bottom_%i.mat',filedirectory, file_n),'file_r')
        end
        x_from_rot = [x_from_rot(:); x_non_ima(:)];
        y_from_rot = [y_from_rot(:); y_non_ima(:)];
    end
    toc
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