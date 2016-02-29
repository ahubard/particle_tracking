%function [fit_x, fit_y, phi_x, phi_y, A_x, A_y] = find_edges(En)
%Find the edges of figures, specifically of the drum

 file_for_center = sprintf('/Users/Aline/Documents/Research/MATLAB/Avalanches/To_get_center_%i.mat',En);
 load(file_for_center);
 nr = length(Total_rotation_steps);

%% Use function edge with 'Canny' method
TOPLEFT = 2*ones(2,nr);
TOPRIGHT = 2*ones(2,nr);
LEFT = zeros(2,nr);
RIGHT = zeros(2,nr);
for ii = 1:nr
    ima = rot_images(:,:,ii);
    ime = edge(ima,'Canny');
    [iy, ix] = find(ime > 0);
    
    %Find top points
    inside = find(ix > 25 & ix < 200);
    ixi= ix(inside);
    iyi = iy(inside);
    auxx = ixi(iyi == 2);
    TOPLEFT(1,ii) = auxx(1);
    
    inside = find(ix > 1050 & ix < 1110);
    ixi= ix(inside);
    iyi = iy(inside);
    auxx = ixi(iyi == 2);
    TOPRIGHT(1,ii) = auxx(end);
    
    %Find lower points
    inside = find(ix > 300 & ix < 400);
    ixi= ix(inside);
    iyi = iy(inside); 
    [sy, ord] = sort(iyi);
    LEFT(1,ii) = ixi(ord(end));
    LEFT(2,ii) = sy(end);
    
    inside = find(ix > 900 & ix < 1000);
    ixi= ix(inside);
    iyi = iy(inside); 
    [sy, ord] = sort(iyi);
    RIGHT(1,ii) = ixi(ord(end));
    RIGHT(2,ii) = sy(end);
%     if (mod(ii,100) == 0)
%         imagesc(ime);
%         hold on;
%         plot(TOPLEFT(1,ii),TOPLEFT(2,ii),'o',TOPRIGHT(1,ii),TOPRIGHT(2,ii),'o',RIGHT(1,ii),RIGHT(2,ii),'o',LEFT(1,ii),LEFT(2,ii),'o')    
%         drawnow
%         pause()
%         hold off
%     end
    
end

%Find center (xot,yot) point that is at same distance of points TOPLEFT,
%TOPRIGHT,RIGHT and LEFT.
%% To use more points for x position. 

xot = (TOPLEFT(1,:) + TOPRIGHT(1,:))/2;

yot1 = (TOPLEFT(2,:)+RIGHT(2,:))/2 + (TOPLEFT(1,:)-RIGHT(1,:)).*...
    (TOPLEFT(1,:)+RIGHT(1,:)-2*xot)./(2*(TOPLEFT(2,:)-RIGHT(2,:)));

yot2 = (TOPRIGHT(2,:)+LEFT(2,:))/2 + (TOPRIGHT(1,:)-LEFT(1,:)).*...
    (TOPRIGHT(1,:)+LEFT(1,:)-2*xot)./(2*(TOPRIGHT(2,:)-LEFT(2,:)));
yot = (yot1+yot2)/2;

% xot = xot(400:end);
% yot = yot(400:end);
% Total_rotation_steps = Total_rotation_steps(400:end);

% Do a first order fourier fit. Values of x vary more so fit is better.
fit_x = fit(Total_rotation_steps',xot','Fourier1');
%% Use the fit_x to exclude the outcast points and get a more accurate
%frequency and amplitude. 
error_fit = xot - fit_x(Total_rotation_steps)';
include = find(abs(error_fit) < 2*std(error_fit));


fit_x = fit(Total_rotation_steps(include)',xot(include)','Fourier1');
fit_y = fit(Total_rotation_steps(include)',yot(include)','Fourier1');


% Put use  xot = xo+A_xsin(wth+phi_x) and yot = yo+A_ycos(wth+phi_y)
phi_x = atan(fit_x.a1/fit_x.b1);
phi_y = atan(-fit_y.b1/fit_y.a1);
A_x = fit_x.a1/sin(phi_x);
A_y = fit_y.a1/cos(phi_y);


%save(file_for_center,'fit_x','fit_y','phi_x','phi_y','A_x','A_y','-append');

