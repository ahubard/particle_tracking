function [bk1 bk2] = getbackground(En,ni,folder,order)
% En is the experiment series number
% ni number of file
% folder number of aline in Poincare
% order indicates if you are calculating bk1 or bk2.
%% Get background goes over images to get the max of each one and normalize by the average light of each image.

%% LOAD IMAGE
fn=sprintf('/aline%i/rotdrum%i/o%02d/onestep%02d_%05d.mat',folder,folder,En,En,ni);
load(fn,'IMA');
bgfile = sprintf('/aline%i/rotdrum%i/o%i/back%i.mat',folder,folder,En,En);




%% Define variables
D = 10;
Ne = size(IMA,3);
bk1 = IMA(:,:,1)*0;
bk2 = bk1;
se = strel('disk',D);
cutoff = 33;  %Taken from histogram

if (order == 2)
    load(bgfile,'bk1');
end
%% Main loop over images
for nframe=1:Ne
    
    im = IMA(:,:,nframe);
    bwim = (imopen(im,se)>cutoff); % bw image. 0 in particles area 1 particle free region
    freezone = bwim.*im;
    norm_coef = sum(freezone(:))/sum(bwim(:));
    im = im/norm_coef;
    
    if (order == 1)
        bk1 = max(bk1,im);
    else
        bk2(im>0)=max(bk2(im>0),bk1(im>0)-im(im>0));
    end
    
    
    
    
    
end

save(fn,'bk1','bk2','-append');


end
