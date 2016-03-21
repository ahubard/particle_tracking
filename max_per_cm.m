yc = 780-ycm(2:end);
ava = avalanche(2:end);
yc = yc(ava > 13)*g*d*m/D*mean(N_particles);
ava = ava(ava > 13)*g*d*m/D;
plot(yc,ava,'.');

