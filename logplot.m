function [xh, yh, expn, normv, p] = logplot(variable,minnorm,N_bin,x_left,x_right)

%creates histogram with logartimic bins

xmin = log10(max(minnorm,min(variable)));
xmax = log10(max(variable)+10*minnorm);
x = logspace(xmin,xmax,N_bin);
h = histogram(variable,x,'Normalization','pdf');
yh = h.Values;
xh = h.BinEdges(2:end);
loglog(xh,yh,'.');
i_f = find(xh > x_left,1,'first');
i_l = find(xh < x_right,1,'last');
xx = xh(i_f:i_l);
yy = yh(i_f:i_l);
fit_h = fit(log(xx(:)),log(yy(:)),'Poly1');
expn = fit_h.p1;
normv = exp(fit_h.p2);
p = loglog(xh,yh,'.',xh,normv*xh.^expn);