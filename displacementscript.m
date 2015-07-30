
function count = displacementscript(folder,En)

avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
%avanofile = sprintf('Avanonestep%i.mat',En);
if(exist(avanofile,'file'))
    load(avanofile);
else 
    save(sprintf('Warning. The file: %s does not exist.mat',avanofile));
    error('Error, avanofile does not exist'); 
end
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
%% Check where avalanches start and end. 
changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);
Nbavalanches = length(initialfileindex);
count = 1;
%% Main loop
for na = 1 :Nbavalanches
 count = displacement(folder,En,na,count,initialfileindex(na),finalfileindex(na),git_version);
end
    
save(avanofile,'count','-append');
