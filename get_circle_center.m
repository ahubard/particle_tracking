function   [xo, yo, R] = get_circle_center(folder,En,file_number,xo_ini,yo_ini,R_ini)
%% Fit circle from particle positions

%% Get Initial range of xo, yo and Radius from firstpositions;

xo_range = fix(xo_ini)+(-10:10);
yo_range = fix(yo_ini)+(-10:10);
range_size = length(yo_range);
R_range = fix(R_ini)+(-2:2);


%% Load file with x,y particle positions
file_for_edge = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',...
        folder,folder,En,En,file_number);
load(file_for_edge);
px = pxs(pxs(:) > 0);
py = pys(pxs(:) > 0);
info.Ny = 400;
info.Nx = 1280;

%% Keep only the external points to reduce the noise. 
%And create image BW_ime = 1 if there is a particle there.
[spy, iy] = sort(py);
spx = px(iy);
ch_i = find(diff(spy) > 0);
BW_ime = zeros(info.Ny,info.Nx);
start_point = 1;

for ii = 1:length(ch_i)
    y = spy(ch_i(ii));
    x_min = min(spx(start_point:ch_i(ii)));
    x_max = max(spx(start_point:ch_i(ii)));
    BW_ime(y,[x_min x_max]) = 1;
    start_point = ch_i(ii)+1;
end

%% Get the best parameters of xo.
%Circle hugh transformation
H = (hugh_circle(BW_ime,xo_range,yo_range,R_range)); %Circle hugh transformation
%Find peaks of the transformation.
[~, I] = max(H(:)); 
[a, ~,~] = ind2sub(size(H),I);
xo = xo_range(a(1)); % Choose main peak for xo. 
rx2 = (px-xo).^2; %Distance of all points to xo

%% Make histograms for the different y values
hist_range = 1:2:R_range(end)+1;
Ri = zeros(length(px),range_size);
di = zeros(range_size,length(hist_range));
ri = zeros(range_size,length(hist_range));

for ii = 1:range_size
    Ri(:,ii) = sqrt(rx2+(py-yo_range(ii)).^2);
    [di(ii,:),ri(ii,:)] = hist(Ri(:,ii),hist_range);
end
%yo is given by the yo_range with more particles at distance 613 from
%origin
[~, iy] = max(di(:,length(hist_range)-1));
%% Keep only external points using (xo, yo) as center and find best fit circle. 
px = px(Ri(:,iy) > R_range(1)-1);
py = py(Ri(:,iy) > R_range(1)-1);
circle_parameters = CircleFitByPratt([px py]);

xo = circle_parameters(1);
yo = circle_parameters(2);
R = circle_parameters(3);



