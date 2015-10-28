function diff_C_Mass = surface_avalanche(folder,En)
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');

%Load needed file information. 
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile,'yo');
fna = sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(fna, 'Avalanche_time','Displacement_File_nb','Number_Avalanches');  
file_CM = sprintf('%sSurface_CM_%i.mat',filedirectory,En);

D =10;
Nb_bins = 123;
diff_CMass_t = cell(1,Number_Avalanches);    
displacement_file = unique(Displacement_File_nb);
number_files = length(displacement_file);
avalanche_counter = 1;

for ii = 1 : number_files
    nf = displacement_file(ii);
    filekernel =sprintf('Displacement_%i',nf);
    fnt =sprintf('%s%s.mat',filedirectory,filekernel);
    load(fnt,'PX','PY');
    t1 = Avalanche_time{1,nf};
    t2 = Avalanche_time{2,nf};
    nb_avalanches = length(t1);
    
    for na = 1:nb_avalanches
        [~, ~, xy_massi, ii_surfacei, bindexi] = ...
            estimate_angle(PX(:,t1(na)),PY(:,t1(na)),D,yo,Nb_bins, 1);
        dummy_aux = zeros (Nb_bins, t2(na)-t1(na)+1);
        for t = t1(na):t2(na)
             dummy_aux(:,t) = compare_C_Mass(PX(ii_surfacei,t),PY(ii_surfacei,t),bindexi,xy_massi,Nb_bins);    
        end
        diff_CMass_t{avalanche_counter} = dummy_aux;
        avalanche_counter = avalanche_counter+1;
        
    end
    
end
        
 save(file_CM,'diff_CMass_t','Avalanche_time','Displacement_File_nb');
       
        
        






    

 
    
    

    
