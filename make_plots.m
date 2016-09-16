savefiles = false;

%% Experimental info
fps = 694.4444;
D = 10;
alpha = 29;
g = 9.8;
d = 1.19e-3;
m = 6.9e-9;
%% Piles that didnt spaned the whole system
%.882
i_stay = find(Pile_Profile_Change <.9 & Num_boundary(1,:)<1.7e4 & Num_boundary(2,:) <25);
i_go = setdiff(1:length(Pile_Profile_Change),i_stay);


%% Energy distribution plots
%Energy distribution Complete plot 
delta_U = [CHANGE_U EXTRA_U];
U = -delta_U*d/D*m;
[xPW, yPW] = logplot(U,1e-12,36,1e-11,1.1e-8,0);
pW = figure(1);
pW.Color = [1 1 1];
pW = loglog(xPW,yPW,'o');

axis([.9*xPW(1) xPW(end-1) min(yPW)*.9 4*max(yPW)]);
axP = gca;
axP.LineWidth = 2;
axP.FontSize = 20;
axP.Box = 'on';
axP.XTick = [1e-12 1e-11 1e-10 1e-9 1e-8 1e-7];
axP.YTick = [1e5 1e6 1e7 1e8 1e9 1e10 1e11];
%labels
xlabel('\DeltaU [J]','FontSize',20);
ylabel ('p(\DeltaU)','FontSize',20);
%Data display
pW(1).MarkerSize = 7;
pW(1).Color = [.5 .7 .5];
pW(1).MarkerFaceColor = [.5 .8 .6];
P_leg = legend('Complete Experimental Data');
P_leg.Box = 'off';
P_leg.FontSize = 16;

%% Distribution of avalanches that are not affected by the system size.

delta_U = [CHANGE_U(i_stay) EXTRA_U];
U = -delta_U*d/D*m;
figure(2)
clf
[xP, yP, expo, normv, p] = logplot(U,1e-12,36,1e-11,1.1e-8,1);
%open figure and plot
pP = figure(2);
pP.Color = [1 1 1];
pP = loglog(xP,yP,'o',xP(3:end-1),normv*xP(3:end-1).^expo);
axis([.9*xP(1) 1.1*xP(end-1) min(yP)*.9 4*max(yP)]);
axP = gca;
axP.LineWidth = 2;
axP.FontSize = 20;
axP.Box = 'on';
axP.XTick = [1e-12 1e-11 1e-10 1e-9 1e-8 1e-7];
axP.YTick = [1e5 1e6 1e7 1e8 1e9 1e10 1e11];
%labels
xlabel('\DeltaU [J]','FontSize',20);
ylabel ('p(\DeltaU)','FontSize',20);
%Data display
pP(1).MarkerSize = 7;
pP(1).Color = [.1 .7 .5];
pP(1).MarkerFaceColor = [.1 .8 .6];
pP(2).LineWidth = 3;
pP(2).Color = [0 0 1];
messa = ('p(\DeltaU) ~ \DeltaU');
P_leg = legend('Experimental Data',sprintf('%s^{%.2f}',messa,expo));
P_leg.Box = 'off';
P_leg.FontSize = 16;


%% Distribution of avalanches affected by the system size.
delta_U = CHANGE_U(i_go) ;
U = -delta_U*d/D*m;
figure(3)
[xPG, yPG] = logplot(U(2:end),1e-12,36,1e-11,1.1e-8,0);
pG = figure(3);
pG.Color = [1 1 1];
pG = loglog(xPG,yPG,'o');

axis([.9*xPW(1) 1.1*xPW(end-1) min(yPW)*.9 4*max(yPW)]);
axP = gca;
axP.LineWidth = 2;
axP.FontSize = 20;
axP.Box = 'on';
axP.XTick = [1e-12 1e-11 1e-10 1e-9 1e-8 1e-7];
axP.YTick = [1e5 1e6 1e7 1e8 1e9 1e10 1e11];
%labels
xlabel('\DeltaU [J]','FontSize',20);
ylabel ('p(\DeltaU)','FontSize',20);
%Data display
pG(1).MarkerSize = 7;
pG(1).Color = [.6 .7 .5];
pG(1).MarkerFaceColor = [.6 .8 .6];
P_leg = legend('Experimental Data');
P_leg.Box = 'off';
P_leg.FontSize = 16;


%% First the log plot of the potential energy.
%Create the variables and fit
figure(1)
clf
[xP, yP] = logplot(-U,1e-11,35,5e-11,8e-7);
[xP_n, yP_n, expn_P, normv] = logplot(-U(ii_non_spaning),1e-11,35,5e-11,8e-7);
%open figure and plot
pP = figure(1);
pP.Color = [1 1 1];
pP = loglog(xP,yP,'.',xP_n,yP_n,'.',xP(3:end-1),normv*xP(3:end-1).^expn_P);
%Data display
pP(1).MarkerSize = 25;
pP(1).Color = [.1 .7 .5];
pP(2).MarkerSize = 25;
pP(2).Color = [.5 0 1];
pP(3).LineWidth = 1.2;
pP(3).Color ='k';
%Axis
axis([xP(1)/2 2*xP(end) min(yP)/2 3*max(yP_n)]);
axP = gca;
axP.FontSize = 20;
axP.Box = 'on';
%labels
xlabel('\DeltaU [J]','FontSize',20);
ylabel ( 'p(\DeltaU)','FontSize',20);
%Legend
messa = ('p(\DeltaU) ~ \DeltaU');
P_leg = legend('Experimental Data', 'Data minus boundary efects',sprintf('%s^{%.2f}',messa,expn_P));
P_leg.Box = 'off';
P_leg.FontSize = 16;



%% Duration of avalanches
%Create the variables and fit

figure(4);
clf
[xT, yT] = logplot(T,1/35,41,1/30,4,0);
[xT_n, yT_n, expn_T, normv] = logplot(T(ii_non_spaning),1/35,41,0.05,4);
%open figure and plot
pT = figure(2);
pT.Color = [1 1 1];
pT = loglog(xT,yT,'.',xT_n,yT_n,'.',xT(3:end-1),normv*xT(3:end-1).^expn_T);
%Data display
pT(1).MarkerSize = 25;
pT(1).Color = [.1 .7 .5];
pT(2).MarkerSize = 25;
pT(2).Color = [.5 0 1];
pT(3).LineWidth = 1.2;
pT(3).Color ='k';
%Axis
axis([xT(1)/2 2*xT(end) min(yT)/2 3*max(yT_n)]);
axT = gca;
axT.FontSize = 20;
axT.Box = 'on';
%labels
xlabel('T [s]','FontSize',20);
ylabel ( 'p(T)','FontSize',20);
%Legend
messa = ('p(T) ~ T');
T_leg = legend('Experimental Data', 'Data minus boundary efects',sprintf('%s^{%.2f}',messa,expn_T));
T_leg.Box = 'off';
T_leg.FontSize = 16;


%% loglog of T vs U
figure(3);
clf
pTU = figure(3);
pTU.Color = [1 1 1];
escaling_f = 5e3;
pTU = loglog(T(ii_non_spaning),-Size_Potential(ii_non_spaning),'.', ...
    T,escaling_f*T.^((-expn_T-1)/(-expn_P-1)));
pTU(1).Color = [.1 .7 .5];
pTU(2).Color = 'k';
pTU(2).LineWidth = 1;

xlabel(' T [s]');
ylabel('\DeltaU [s]');
ax = gca;
ax.FontSize = 20; 
messag = '\DeltaU ~ T'; 
leg = legend(['Experimental Data corrected ' char(10)  ' for boundary effects'],...
    sprintf('%s^{%.2f}',messag,((-expn_T-1)/(-expn_P-1))));
leg.Location = 'northwest';
leg.Box = 'off';

%% Average shape
it = setdiff(1:16,[7 9 13]);
figure(4)
clf
clf
AS = figure(4);
AS.Color = [1 1 1];
AS = plot(average_shape(:,it)/1e-11,'-','LineWidth',1);
axis([1 101 0 5]);
ax = gca;
ax.XTick = 1:20:101;
ax.XTickLabel = 0:.2:1;
ax.YTick = (1:1:5);
ax.FontSize = 20;
xlabel ('t/T');
ylabel('<\DeltaU_{t}>_{T}');
LEG = regexp(sprintf('T = %.2f s#', duration), '#', 'split');
LEG(end) = [];
LEG = legend(LEG);
LEG.Location = 'eastoutside';

%% Collapsed Average Shape


xmin = log10(.08);
xmax = log10(1);
x = logspace(xmin,xmax,21);
figure(5)

clf
ASN = figure(5);
ASN.Color = [1 1 1];
ASN = plot(average_shape(:,it)./repmat(x(it),101,1)/1e-10,'LineWidth',1);
hold on;
%ASN(14) = plot(mean(average_shape(:,it)./repmat(x(it),101,1)/1e-10,2),'k','LineWidth',1.5);
axis([1 101 0 1]);
ax = gca;
ax.XTick = 1:20:101;
ax.XTickLabel = (0:.2:1);
%ax.YTick = (.2:.2:1)*1e-10;
ax.FontSize = 20;
xlabel ('t/T');
ylabel('T^{1-\sigma}<\DeltaU_{t}>_{T}');
LEG = regexp(sprintf('T = %.2f s#', duration), '#', 'split');
LEG(end) = [];
LEG = legend(LEG);
LEG.Location = 'eastoutside';
hold off

%% Save the images
if(savefiles)
file_name1 = '/Users/Aline/Documents/Research/Presentations/APS_2016/Dist_Pot.pdf';
save_fig(1,file_name1);
file_name2 = '/Users/Aline/Documents/Research/Presentations/APS_2016/Dist_Duration.pdf';
save_fig(2,file_name2);
file_name3 = '/Users/Aline/Documents/Research/Presentations/APS_2016/EnergyvsDuration.pdf';
save_fig(3,file_name3);
file_name4 = '/Users/Aline/Documents/Research/Presentations/APS_2016/AVESHAPE.pdf';
save_fig(4,file_name4);
file_name4 = '/Users/Aline/Documents/Research/Presentations/APS_2016/AVESHAPEcoll.pdf';
save_fig(5,file_name4);
file_name10 = '/Users/Aline/Documents/Research/Presentations/APS_2016/Dist_Potsans.pdf';
save_fig(10,file_name10);
end