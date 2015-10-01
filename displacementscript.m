
function count = displacementscript(folder,En)

avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
%avanofile = sprintf('Avanonestep%i.mat',En);
if(exist(avanofile,'file'))
    load(avanofile);
else 
    save(sprintf('Warning. The file Avanonestep%i does not exist.mat',En));
    error('Error, avanofile does not exist'); 
end
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
alpha = 29; %Angle of the camera from experiment
maxT = 3248; %got maximal duration from previous data.  
%% Check where avalanches start and end. 
changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);
Nbavalanches = length(initialfileindex);
count = 1;
%% Main loop
for na = 1 :Nbavalanches
 count = displacement(folder,En,na,count,initialfileindex(na),finalfileindex(na),alpha,git_version);
end
    
save(avanofile,'count','alpha','maxT','-append');
