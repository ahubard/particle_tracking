folder = 1;
En = 104;
%avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
avanofile = sprintf('Avanonestep%i.mat',En);
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
%% Check where avalanches start and end. 
changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);
%Nbavalanches = length(initialfileindex);
Nbavalanches = 10;
%% Main loop
for na = 1 :Nbavalanches
 displacement(folder,En,na,initialfileindex(na),finalfileindex(na),git_version);
end
    