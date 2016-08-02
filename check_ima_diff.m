function  [ima_diff, same_rot, file_same_rot1, file_same_rot2] = check_ima_diff(folder,En)

%% Get rotations info
filedirectory = sprintf('/aline%i/rotdrum%i/o%02d/',folder,folder,En);
avanofile = sprintf('%sAvanonestep%i.mat',filedirectory,En);
load(avanofile,'avan');
n_files = size(avan,2);
file_exist = zeros(1,n_files);
with_undersquare = zeros(1,n_files);
%% Loop over the n_files to see how many image files are available
for ii = 1:n_files
    imagefile1 = sprintf('%sonestep%i%05i.mat',filedirectory,En,ii);
    if (exist(imagefile1,'file'))
        file_exist(ii) = 1;
    else
        imagefile2 = sprintf('%sonestep%i_%05i.mat',filedirectory,En,ii);
        
        if(exist(imagefile2,'file'))
            file_exist(ii) = 1;
            with_undersquare(ii) = 1;
            
        end
    end
end

 %% decide wich one is the right file name
 if(sum(with_undersquare) == 0)
     filebase = sprintf('%sonestep%i',filedirectory,En);
 else
     filebase = sprintf('%sonestep%i_',filedirectory,En);
 end
 
%% Get the rotation step of all the existing files and find when they were
% recorded after the same rotation.
file_num = find(file_exist);

same_rot = find(diff(avan(2,file_num))==0);
file_same_rot1 = file_num(same_rot);
file_same_rot2 = file_num(same_rot+1);
loop_diff = avan(1,file_same_rot2)-avan(1,file_same_rot1); %If one is continuos;
ns_rot = length(file_same_rot1);
ima_diff = zeros(1,ns_rot);
%% Main loop to get ima_diff
for nf = 1:ns_rot
    fn1 = file_same_rot1(nf);
    fn2 = file_same_rot2(nf);
    I1 = load(sprintf('%s%05i.mat',filebase,fn1),'IMA');
    I2 = load(sprintf('%s%05i.mat',filebase,fn2),'IMA');
    ima_diff(nf) = sum(sum((I1.IMA(:,:,351)-I2.IMA(:,:,1)).^2));
end

save (sprintf('%sIMA_DIFF_%i.mat',filedirectory,En));   

       
       
