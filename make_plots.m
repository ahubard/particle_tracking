%% Script to make the plots for the avalanches. 

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
P_leg = legend('Experimental Data', 'Data minus boundary efects',sprintf('%s^{%.1f}',messa,expn_P));
P_leg.Box = 'off';
P_leg.FontSize = 16;

file_name = '/Users/Aline/Documents/Research/Presentations/APS_2016/Dist_Pot.pdf';
save_fig(1,file_name);

%% Duration of avalanches
%Create the variables and fit

figure(2);
clf
[xT, yT] = logplot(T,1/35,41,1/30,4);
[xT_n, yT_n, expn_T, normv] = logplot(T(ii_non_spaning),1/35,41,0.04,3);
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
T_leg = legend('Experimental Data', 'Data minus boundary efects',sprintf('%s^{%.1f}',messa,expn_T));
T_leg.Box = 'off';
T_leg.FontSize = 16;
file_name = '/Users/Aline/Documents/Research/Presentations/APS_2016/Dist_Duration.pdf';
save_fig(2,file_name);

%% 
