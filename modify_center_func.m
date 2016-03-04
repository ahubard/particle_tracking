%% Get rid of outliers and use that the frequency is already given for the fit. 

[git_version, ~] = evalc('system(''git describe --dirty --alway'')');


%% Use images and angle of rotation to find how the center of rotation moves
%filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
%file_for_center = sprintf('%sCenter_%i.mat',filedirectory,En);
%file_for_center = '/Users/Aline/Documents/Research/MATLAB/Avalanches/
%file_for_center
%load(file_for_center,'git_version','nb_rotation_steps','xot','yot','Rt');
w = 0.0032;
%th = Total_rotation_steps*w;
th = nb_rotation_steps*w;
%% First get rid of outliers by doing a linear aproximattion
a0_x = mean(xot);
include_x = find(abs(xot-a0_x) < 6);
a0_y = mean(yot);
include_y = find(abs(yot-a0_y) < 6);


ft = fittype({'sin(x)','cos(x)','1'},'coefficients',{'b','c','a0'});
ft_x = fit(th(include_x)',xot(include_x)',ft);
% get xot = Asin(th+phi_x)+a0
phi_x = atan(ft_x.c/ft_x.b);
a0_x = ft_x.a0;
A_x = ft_x.b/cos(phi_x);

%% Used first to get the relation btw xot and yot
ft_y = fit(th(include_y)',yot(include_y)',ft);
%get xot = Asin(th+phi_x)+a0
phi_y = atan(-ft_y.b/ft_y.c);
a0_y = ft_y.a0;
A_y = ft_y.c/cos(phi_y);
%% Use info from previuos fit for files En>100
Amplitud_coef = A_y/A_x;
delta_phi = phi_y - phi_x;
