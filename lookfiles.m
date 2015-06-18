En =103;
folder  =1;
avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
load(avanofile,'avan','navfile');
changefileindex = find(diff(avan(1,navfile))>1);
initialfileindex = navfile([1 changefileindex(1:end-1)+1]);
finalfileindex = navfile(changefileindex);


for ii= 1:length(initialfileindex)
tf =sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,ii);
naf = sprintf('/aline%i/rotdrum%i/o%02d/NoAvalanche_%i.mat',folder,folder,En,ii);
trackf(ii) =exist(tf,'file');
noava(ii) =exist(naf,'file');
end

nofile = find(trackf+noava == 0);
firstfile = initialfileindex(nofile);
finalfile = finalfileindex(nofile);
length(firstfile)


for ii = 1:length(firstfile);
    NEn = nofile(ii);
     j(ii)=batch(sched,'mytrack',1,{folder,En,NEn,firstfile(ii),finalfile(ii),D},'Filedependencies',{'stickfiles.m','adjacent.m','assignmentoptimal.m'});
end

    
    