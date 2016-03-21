function save_fig(n,file_name)
figure(n);
ax = gca;
ax.FontSize = 20;
ax.Box = 'on';
%Fit fig to space
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1)-ti(3);
ax_height = outerpos(4) - ti(2)-ti(4);
ax.Position = [left bottom ax_width/1.01 ax_height/1.01];
print(file_name,'-dpdf');
% fig = gcf;'-dpdf'
% saveas(fig,file_name);
