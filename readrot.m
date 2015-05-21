function [im, HEAD, TAIL]=readrot(fn,nx,ny,N,nhead,ntail)

%% File to read the frames from a the bin files of the avalanches
% fn is the name of the file
%nx, ny size of desired image
% frame number to read
% nhead and ntail, head and tail of each image
%% Open the file
% if (N==1)
%     fprint('corrupted image, start form the second')
% end

fid=fopen(fn);     % Handle to file/
%% Set position and size to read to read; 
framesize=nhead+ntail+nx*ny;       % Frame size with head and tail
offset=(N-1)*(nhead+ntail+nx*ny);  %From whic byte start reading
fseek(fid,offset,'bof');            %Move to position offset in bytes in the bin file with handle fid
raw=fread(fid,framesize);           %Read one frame
im=(reshape(raw(nhead+1:end-ntail),nx,ny))';   %Reshape to a nx ny matrix
HEAD=raw(1:nhead);
TAIL=raw(end-ntail:end);
fclose(fid);