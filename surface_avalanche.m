function diff_C_Mass = surface_avalanche(folder,En)
%\
avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
%avanofile = sprintf('Avanonestep%i.mat',En);
if(exist(avanofile,'file'))
    load(avanofile);
else 
    save(sprintf('Warning. The file Avanonestep%i does not exist.mat',En));
    error('Error, avanofile does not exist'); 
end
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');

%% Check where avalanches start and end. 
changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);
Nbavalanches = length(initialfileindex);
count = 1;
%% Main loop
for na = 1 :Nbavalanches
 
    filekernel =sprintf('Displacement_%i',na);

    fnt =sprintf('%s%s.mat',filedirectory,filekernel);
    load(fnt,'PX','PY');

    %use avalanche time from avalanche_size files
    
    
    
end
    