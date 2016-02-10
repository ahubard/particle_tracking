function  [average_displacement, total_overlap, XI2]  = compare_rotation(alpha_1)


folder = 1; 
En = 103;


filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);
rotation_center_file = sprintf('%srotation_center_%i.mat',filedirectory,En);
load(rotation_center_file);

angle_aperture = 2.35;
i_final = find(diff(Rotation_step(1,:))~=0);
i_initial = find(diff(Rotation_step(1,:))~=0)+1;d_R = Rotation_step(1,i_final)-Rotation_step(1,i_initial);
nb_rotation_steps_btw_avalanches = Rotation_step(1,i_initial)-Rotation_step(1,i_final);

for ii = 1:length(i_initial)
    before_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,Fn_imafile(i_final(ii)));
    after_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,In_imafile(i_initial(ii)));
    load(before_fn,'pxs','pys','Npf');
    pxb(1:Npf(351),ii) = pxs(1:Npf(351),351);%before rotation, after the avalanche finishes.
    pyb(1:Npf(351),ii) = pys(1:Npf(351),351);
    load(after_fn,'pxs','pys','Npf')
    pxa(1:Npf(1),ii) = pxs(1:Npf(1),1);%after rotation, before the avalanche starts.
    pya(1:Npf(1),ii) = pys(1:Npf(1),1);
end


 


%% Check to see if it matches.
load(rotation_center_file);
%% Find rotation angle
%alpha_1 = .00315; %initial angle. 
alpha = nb_rotation_steps_btw_avalanches*alpha_1;
alpha = -repmat(alpha,size(pxb,1),1);
tran_pxb = (pxb-xo).*cos(alpha)-(pyb-yo).*sin(alpha)+xo;
tran_pyb = (pxb-xo).*sin(alpha)+(pyb-yo).*cos(alpha)+yo;


X_non_rot = zeros(2210,size(pxa,2));
Y_non_rot = zeros(2210,size(pxa,2));
X_rot = zeros(2210,size(pxa,2));
Y_rot = zeros(2210,size(pxa,2));
overlaping_number = zeros(1,size(pxa,2));
%% Find overlap points
for ii = 1:size(pxa,2)
[~, trivialbondt1,trivialbondt2] = adjacent(pxa(:,ii),pya(:,ii),tran_pxb(:,ii),tran_pyb(:,ii),2);

overlaping_number(ii) = length(trivialbondt1); 

X_non_rot(1:overlaping_number(ii),ii) = pxa(trivialbondt1,ii);
X_rot(1:overlaping_number(ii),ii) = tran_pxb(trivialbondt2,ii);

Y_non_rot(1:overlaping_number(ii),ii) = pya(trivialbondt1,ii);
Y_rot(1:overlaping_number(ii),ii) = tran_pyb(trivialbondt2,ii);
end


total_diff = sum((X_non_rot-X_rot).^2+(Y_non_rot-Y_rot).^2);
total_overlap = sum(overlaping_number);
XI2 = sum(total_diff);
average_displacement = XI2/total_overlap;
%y1y2 = [y1 y2];

% 
% [low_x il_x] = min(x1x2);
% vx(1) = low_x + D*4;
% vy(1) = y1y2(il_x) +D*4;
% 
% [low_y il_y] = min(y1y2);
% vy(2) = low_y + D*10;
% vx(2) = x1x2(il_y);
% 
% [max_x im_x] = max(x1x2);
% vx(3) = max_x - D*4;
% vy(3) = y1y2(im_x) +D*4;
% 
% [max_y im_y] = max(y1y2);
% vy(4) = max_y;
% vx(4) = x1x2(im_y);
% 
% low_y = -400;
% max_y = yo+R;
% ii_b = find(tran_pxb(:,ii) < max_x  &  tran_pxb(:,ii) > low_x  & ...
%     tran_pyb(:,ii) < max_y & tran_pyb(:,ii) > low_y);
% ii_a = find(pxa(:,ii) < max_x  &  pxa(:,ii) > low_x  &  pya(:,ii) < max_y &...
%     pya(:,ii) > low_y);
%  [~, trivialbondt1,trivialbondt2] = adjacent(pxa(ii_a,ii),pya(ii_a,ii),...
%      tran_pxb(ii_b,ii),tran_pyb(ii_b,ii),2);

% X_non_rot(1:overlaping_number,ii) = pxa(ii_a(trivialbondt1),ii)';
% Y_non_rot(1:overlaping_number,ii) = pya(ii_a(trivialbondt1),ii)';
% X_rot(1:overlaping_number,ii)  = tran_pxb(ii_a(trivialbondt2),ii)';
% Y_rot(1:overlaping_number,ii)  = tran_pyb(ii_a(trivialbondt2),ii)';
end
