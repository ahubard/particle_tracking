function avaseq = makeavaseq(folder,En)
avanofile = sprintf('/aline%i/rotdrum%i/o%i/Avanonestep%i.mat',folder,folder,En,En);
load(avanofile);

%% Find maximum duration of avalanche
nbT = length(nbtracked);
filesize = zeros(1,nbT);
first = zeros(1,nbT);
last = zeros(1,nbT);
for nb = 1:nbT
     fnn =sprintf('/aline%i/rotdrum%i/o%02d/Trackeds_%i.mat',folder,folder,En,nbtracked(nb));
     load(fnn,'initial','final');
     first(nb) = initial;
     last(nb) = final;
end



%% Load files to create seqofavalanches 
 avaseq = zeros(1,maxavaduration*nr);  %nr is the number of rotations in the experiment.
 
 