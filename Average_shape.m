function [nbavalanchesforshape, averageshape] = Average_shape(T1,T2,increment)
%average shape of avalanches of size T T1<T<T2
%% Load needed files
basename = sprintf('\\\\POINCARE\\aline_share\\rotdrum1\\');  %when running local

averageshape = zeros(1,101);
nbavalanchesforshape = 0;

for exprun = 14:20
    fname = sprintf('o%i\\AVALANCHE_%i.mat',exprun,exprun);
    filename = strcat(basename,fname);
    load(filename);
    iava = find ((avalancheduration<T2) .* (avalancheduration>T1));
    %nbavalanchesforshape = nbavalanchesforshape + length(iava);
    
    for nava = 1:length(iava)
        na = iava(nava);
        avalanche = Kinetic_energy(avainit(na):avaendt(na),avaini(na));
        deltat = avaendt(na)-avainit(na);
        interpolateava = interp1(([0 increment/2+(0:deltat) deltat+increment]) /(deltat+increment),[0 avalanche' 0],(0:.01:1),'pchip');  
        averageshape = averageshape + interpolateava/max(interpolateava);
        nbavalanchesforshape = nbavalanchesforshape + 1;
        plot(0:.01:1,averageshape/nbavalanchesforshape);
        drawnow
    end
end

plot(0:.01:1,averageshape/max(averageshape))