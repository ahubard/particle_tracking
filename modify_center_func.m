function [Amp a0_x a0_y phi_x] = modify_center_func(folder,En)
%% Get rid of outliers and use that the frequency is already given for the fit. 
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');


%% Use images and angle of rotation to find how the center of rotation moves
filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
%filedirectory = sprintf('/Users/Aline/Documents/Research/MATLAB/particle_tracking/');

file_for_center = sprintf('%sCenter_%i.mat',filedirectory,En);


load(file_for_center,'git_version','nb_rotation_steps','xot','yot','Rt');
w = 0.0032;
%th = Total_rotation_steps*w;
th = nb_rotation_steps*w;
%% First get rid of outliers by doing a linear aproximattion
a0_x = mean(xot);
include_x = find(abs(xot-a0_x) < 6);
% a0_y = mean(yot);
% include_y = find(abs(yot-a0_y) < 6);


% ft = fittype({'sin(x)','cos(x)','1'},'coefficients',{'b','c','a0'});
% ft_x = fit(th(include_x)',xot(include_x)',ft);
% % get xot = Asin(th+phi_x)+a0
% phi_x = atan(ft_x.c/ft_x.b);
% a0_x = ft_x.a0;
% A_x = ft_x.b/cos(phi_x);

%% Used first to get the relation btw xot and yot
% ft_y = fit(th(include_y)',yot(include_y)',ft);
% %get xot = Asin(th+phi_x)+a0
% phi_y = atan(-ft_y.b/ft_y.c);
% a0_y = ft_y.a0;
% A_y = ft_y.c/cos(phi_y);
% Amplitud_coef = A_y/A_x;
% delta_phi = phi_y - phi_x;
%% Use that all files have Ay = Ax and phi_y = phi_x + .2;

%% Use info from previuos fit for files En < 100
if(En <100)
    Amp = 1.15;
    a0_x = 633.5;
    a0_y = 203;
    R = 579.9;
else
    Amp = 1.76;
    a0_x = 643.7;
    a0_y = 159.4;
    R = 617.7;
    
end
    
%% Rescale data for new fit and find the phase change
n_xot = (xot - a0_x)/Amp;
n_yot = (yot - a0_y)/Amp;
ft = fittype({'sin(x)','cos(x)'},'coefficients',{'b','c'});
n_fit = fit(th(include_x)',n_xot(include_x)',ft);
b = n_fit.b; %for correction of phase if amplitude is positive
phi_x = atan(n_fit.c/n_fit.b) + pi*(sign(b)-1)/2;
phi_x = mod(phi_x,2*pi);
dphi = 0.1;
% figure(1);
% plot(th(include_x),n_xot(include_x),'.',th,sin(th+phi_x));
% drawnow;
% figure(2);
% plot(th(include_x),n_yot(include_x),'.',th,cos(th+phi_x+dphi));
% drawnow;
% n_yot = (yot - a0_y)/Amp;
% n_th = th + phi_x;
% n_fit = fit(th(include_y)',n_yot(include_y)',ft);
% b_y = n_fit.b;
% d_phi = atan(-n_fit.b/n_fit.c);
%% Save to file
save(file_for_center,'git_version','Amp','a0_x','a0_y','phi_x','w','dphi','-append');
