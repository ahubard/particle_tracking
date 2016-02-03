for En= 21:23
folder = 2; 
fprintf ('Your values are En = %i and folder = %i\n',En,folder);
%fprintf('Press enter to continue\n');
%pause();
avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
load(avanofile);
changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);

%% Launch mytrack to conect particle positions.
for ii = 1:length(initialfileindex)
   j(ii)=batch(sched,'mytrack',1,{folder,En,ii,initialfileindex(ii),finalfileindex(ii),D,mk},'Filedependencies',{'stickfiles.m','adjacent.m','assignmentoptimal.m'});
     
end
end