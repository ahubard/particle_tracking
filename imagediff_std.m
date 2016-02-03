file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);
small_change = find(abs(Avalanche_potential/max(NoParticles_moved))<.05);
low_part = find(Participation(small_change) < 0.1);
tiny_ava = small_change(low_part);
imagefile = zeros(2,length(tiny_ava));

for ni = 1:length(tiny_ava)
in =  Displacement_File_nb(tiny_ava(ni));
fnn = sprintf('/aline%i/rotdrum%i/o%02d/Displacement_%i.mat',...
    folder,folder,En,in);
load(fnn,'initial','final');
imagefile(:,ni)= [initial final];
end

recheck = find((imagefile(2,:)-imagefile(1,:))<2);
ii_imagesaux = imagefile(:,recheck);
[ii_images index] = unique(ii_imagesaux(:));
[fi, ni] = ind2sub(size(ii_imagesaux),index);

imagestd = zeros(1,length(ii_images));

for ii = 1:length(ii_images)
    fno = sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',...
        folder,folder,En,En,ii_images(ii));
    load(fno,'IMA');
    ima_diff = (IMA(:,:,1)-IMA(:,:,351));
    imagestd(ii) = std(ima_diff(:));
end
no_avalanche = find(imagestd > 6); 

[fi ni] = ind2sub(size(ii_imagesaux),index(no_avalanche));

    