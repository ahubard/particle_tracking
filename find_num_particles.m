function [mean_nb_particles, std_nb_particles] = find_num_particles(folder,En,initial,final)

filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
nfiles = final-initial+1;
NPF = zeros(1,nfiles*350);

for ii = initial:final
    fn = sprintf('%spositions%02d_%05d.mat',filedirectory,En,ii);
    load(fn,'Npf');
     NPF(((ii - initial)*350)+1:(ii-initial+1)*350) = Npf(1:350);
end

mean_nb_particles = mean(NPF);
std_nb_particles = std(NPF);

