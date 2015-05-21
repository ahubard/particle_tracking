%% avasnap fermi
%% Controls camera to take continuos pictures and save them when there is an avalanche. 
% while telling the motor to rotate when nothing is happened in certain
% time

%% Initialize motor of drum
p3 = serial('COM3','BaudRate',9600);   %Define port
p3.Terminator='CR/LF';

%open port
   fopen(p3); 

   %Define step and velocity
fprintf(p3,'D8'); fscanf(p3)
fprintf(p3,'I20'); fscanf(p3)
fprintf(p3,'V20'); fscanf(p3)



%% General info

%En=15;   %Experiment number
NP=4000;  %Number of saved files. 
NR=3000;
NS=350;    %Number of frames to save/check at a time
NM= 12250;   %Number of frames to pass with no activity before rotating.
ssave=10;   %number of rotations to let go with out checking or saving

% %%

% NP=10;  %Number of saved files. 
% NR=30;
% NS=350;    %Number of frames to save/check at a time
% NM=1200;   %Number of frames to pass with no activity before rotating.
% ssave=1;   %number of rotations to let go with out checking or saving
%% Camera info 

skip=1; %indicates the camera not to skips frames to save in the buffer

% Load camera libary
libname='XCLIBW64';
hlibname='XCLIBW64_IMPORT.H';

if not(libisloaded(libname))
  loadlibrary(libname,hlibname);
end

% open camera with specifications in file imseq
calllib(libname,'pxd_PIXCIopen','', '', 'imseq4.fmt');
% Grab the specifications of pictures
info.Nx=calllib(libname,'pxd_imageXdim');  %Length of picture in pixels
info.Ny=calllib(libname,'pxd_imageYdim');  %Height in pixels
info.Nclr=calllib(libname,'pxd_imageCdim');  %color 1 for monochorme
info.Nbit=calllib(libname,'pxd_imageBdim'); %Number of bits per pixel
info.Nbuf=calllib(libname,'pxd_imageZdim');  %Total amount of pictures than can be saved in the buffer. 
info.fr=calllib(libname,'pxd_SV2112_getCtrlFrameRate',1);
%close the camera
calllib(libname,'pxd_PIXCIclose');  


 nbuf=floor(info.Nbuf/(NS+1))*(NS+1) %Just if smaller than info.Nbuf
 if (nbuf>info.Nbuf)
     nbuf=info.Nbuf
 end
info.nbuf=nbuf;

 pause(.1)
%% INFO for saved info
avafile=sprintf('Avanonestep%02d.mat',En);
ofn=sprintf('onestep%02d%05d.bin',En,1);  %Image files
finfo=sprintf('infof%02d%05d.mat',En,1);   %info file.
save(finfo,'info');
%savefile1=sprintf('savefile1.mat');   %File to indicate to the other terminal if it should save raw data. 
%savefile2=sprintf('savefile2.mat');   %File to indicate to the other terminal if it should save raw data. 
avan=zeros(3,NR*2);                 % To save in which rotation is the file being saved.
checking=zeros(1,NR*20);            %Check where the camera is after each turn of the loop
%imadif=zeros(1,NR*20);
ims=zeros(info.Nx,info.Ny,NS+1,'uint8'); %Save the frames. 
%% Initialize counters
nr=0;   %Rotation step;
na=0;   %avalanche step
ns=1;   %loop step.

%% Avalanche info
%images at the begging and end of each Ns cycle. 
sizeim=info.Nx*info.Ny;             %image size
im=zeros(info.Nx,info.Ny,'uint8'); %allocate for images
imi=zeros(info.Nx,info.Ny,'uint8'); %allocate first image
imf=zeros(info.Nx,info.Ny,'uint8'); %allocate last image
%cutoff=2.5e6; With out using acerage of image to normalize
cutoff=2e3;  %Divindinge by the average light per pixel of each image
%% Info for motor
 rot=sprintf('+1');  %rotation step size
 move=NM+1;           %Counter of frames with out activity. 
 %motorfile=sprintf('rotation.mat');  %When the file is created the motor moves.

 %% Start camera
 calllib(libname,'pxd_PIXCIopen','', '', 'imseq4.fmt'); %open the camera
%  fprintf('Open camera in another window');
  pause(1);
 calllib(libname,'pxd_goLiveSeq',1,1,nbuf,1,0,skip);  %Load the camera buffers cyclic. 
 %fr =current frame on buffer   % nas number of buffer tha have gone by
 %since the camera started. 
 pause(1.2);
 nac=calllib(libname,'pxd_videoFieldCount',1); frini=calllib(libname,'pxd_capturedBuffer',1) 
 NAC1=nac; %reference or where it started saving pics. 
 fr=NS*3+1
 m1=nac+fr-frini %initial frame of the avalanche count check.
 r2=m1; %initial frame of the rotating count frame. 
 info.nac=m1;


%% Main run while the number of rotation steps is smaller than NS
 while(na<=NP)
     
     %check to see if the motor is ready to move
%      wait=exist('rotation.mat','file');
%      while(wait>0)
%          wait=exist('rotation.mat','file');
%      end
     
%Check to see if the motor moves. 
    if(move>=NM)
        %save(motorfile,'rot');  %save file to tell motor to move
        fprintf(p3,rot);fscanf(p3);
        nr=nr+1     %increment steps
        r1=r2;        %actualize rot counter
       % m1=m2;        %actualize check counter  
    end
    
    %Check last loaded buffer
    m2=calllib(libname,'pxd_videoFieldCount',1);

    %Wait for NS buffers to load
    while((m2-m1)<NS)
        pause(.01);
        m2=calllib(libname,'pxd_videoFieldCount',1);
    end
     
    m1=m1+NS;
    r2=r2+NS;
    %Load  last image
    [ch1 imi]=calllib(libname,'pxd_readuchar',1, mod(fr-1,nbuf)+1, 0, 0, info.Nx, info.Ny, imi, sizeim, 'Gray'); %load first image
    [ch2 imf]=calllib(libname,'pxd_readuchar',1, mod(fr+NS-1,nbuf)+1, 0, 0, info.Nx, info.Ny, imf, sizeim, 'Gray');
    imadif=sum((double(imi(:))/mean(double(imi(:)))-double(imf(:))/mean(double(imf(:)))).^2);
    %imadif=sum((double(imi(:))-double(imf(:))).^2);
    %Check for avalanche
    if(nr>ssave)
    if(imadif>cutoff)
        for nfr=1:NS+1
               [ch2 im]=calllib(libname,'pxd_readuchar',1, mod(fr+nfr-2,nbuf)+1, 0, 0, info.Nx, info.Ny, im, sizeim, 'Gray');
               ims(:,:,nfr)=im;
               im=zeros(info.Nx,info.Ny,'uint8');
        end
        fid=fopen(ofn,'w');
        fwrite(fid,ims,'uint8');
        fclose(fid);
        ims=zeros(info.Nx,info.Ny,NS+1,'uint8');
 

        na=na+1;
        avan(:,na)=[ns; nr; fr];
        r1=r2;
        ofn=sprintf('onestep%02d%05d.bin',En,na+1);  %Image files
        %finfo=sprintf('infof%05d.mat',na+1);   %info file.
    end
    end
    
    
    move=r2-r1;
    checking(ns)=calllib(libname,'pxd_videoFieldCount',1);
    ns=ns+1;
    %info.nac=nac+NS;
    fr=fr+NS;
    %IM(:,:,ns)=double(imi);
    %imi=imf;
    if (mod(na,100)==0)
        save(avafile,'NAC1','avan','checking','na','NS','ns','nr');
    end
    
 end
    
% Nc=length(sprintf('%d',NF));
% fmt=sprintf('[%%0%dd/%%0%dd]',Nc,Nc);
% fprintf('Acquiring...');
% fprintf(1,fmt,0,NF);

save(avafile,'NAC1','avan','checking','na','NS','ns','nr');
fprintf(1,'Done.\n');

calllib(libname,'pxd_goUnLive',1);
%close camera
calllib(libname,'pxd_PIXCIclose');

 

fclose(p3);
%close motor
% rot=sprintf('0');  %rotation step size
% finish=1;
% save(motorfile,'rot','finish');

       %Nf=nbuf-mod(fr-1,nbuf)+1;
%         if(mod(na,3)==0)
%             finishsaving1=exist(savefile1,'file');
%             while(finishsaving1>0)
%                 pause(.01);
%                 finishsaving1=exist(savefile1,'file');
%             end
%             save(savefile1,'ofn','ofnb','fr','NS','Nf');
%         elseif(mod(na,3)==1)
%             finishsaving2=exist(savefile2,'file');
%             while(finishsaving2>0)
%                 pause(.01);
%                 finishsaving2=exist(savefile2,'file');
%             end
%             save(savefile2,'ofn','ofnb','fr','NS','Nf');
%             
%         else
            
   %if(sum((double(imi(:))-double(imf(:))).^2)>cutoff)
       % fprintf(1,'Saving...');
        %check where in the bufffer are you;
%         if (mod(fr-1,nbuf)+1+NS<=nbuf)
        %Save from frame fr to frame fr+NS
       % calllib(libname,'pxd_saveRawBuffers',1,ofn,mod(fr-1,nbuf)+1,mod(fr+NS-1,nbuf),[],0,0,1);
       
        %info.Nf=NS;
%         else
%             %save two files, the one with the last buffers and the one with
%             %the first ones.
%             
%             calllib(libname,'pxd_saveRawBuffers',1,ofn,mod(fr-1,nbuf)+1,nbuf,[],0,0,1);
%             calllib(libname,'pxd_saveRawBuffers',1,ofnb,1,Nf,[],0,0,1);
%         end
%        end

% save(savefile1,'ofn','ofnb','fr','NS','finish','Nf');
% save(savefile2,'ofn','ofnb','fr','NS','finish','Nf');
backup=sprintf('backup%02d.mat',En);
save(backup,'NAC1','avan','checking','na','NS','ns','nr');