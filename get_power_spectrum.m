%Find the power spectrum of each avalanche and get the average. 

Fs = 694444444/1000000; %Sampling frequency, from the camera.
T = 1/Fs;               %Time between frames.

filenumbers = [ 15 16 17 18 19 20 21 22 23 103 104 105 106 107 108 109 ]; %Files that contain the info
Nofiles = length(filenumbers);
kind = 0;
N = 8192;
power_spectrum = zeros(1,N/2);
N_ava = 3392; %from previous runs

%% Get the spectrum for the files of each run
for ii = 1:Nofiles;
    En = filenumbers(ii);
    filename = sprintf('Avalanches_%i_%i.mat',En,kind);
    load(filename,'mat_potential');
    mat_potential = mat_potential';
    L = size(mat_potential,2);      %length of the avalanches in frames.
    t = (0:L-1)*T;                  % time  of the avalanche in sec
    n = 2^nextpow2(L);              %extend length to be a power of two
    N(ii) = n;
    dim = 2;                        %indicates mat_potential is a matrix  
    Y = fft(mat_potential,n,dim);   %Fast fourier transform
    %swicht the transform so it makes sense;
    P2 = abs(Y/n);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    power_spectrum(1:n/2) = power_spectrum(1:n/2) + sum(P1(:,1:n/2))/N_ava;
end
    %loglog(0:(Fs/N):(Fs/2-Fs/N),power_spectrum);
