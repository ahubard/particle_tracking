%% Calls the function avalanche_size, for the different experiments in all 3 kinds. 
%uses a parallel loop to use the cluster and try to reduce the time. 

%% Create the array of different initial positions.
Experiment = [ 15 16 17 18 19 20 21 22 23 103 104 105 106 107 108 109 ];
FOLDER = [1 1 1 1 1 1 2 2 2 1 1 2 2 2 2 2];
LE = length(Experiment);
%% Open matlab in parlel mode
matlabpool open
%regular loop over kind of avalanche separation
for kind = 0:2 
    parfor ie = 1:LE
        avalanche_size(FOLDER(ie),Experiment(ie),kind);
    end
    sprintf('Done with kind = %i',kind);
end
matlabpool close
    