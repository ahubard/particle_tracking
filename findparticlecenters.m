function  keep = findparticlecenters(En,ni,folder,D,w,Cutoff,MinSep)
%% findinrot  find the particles in the rotating drum.
[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
%% Check if it was previosly tracked

track = 1;

%% Check if there is really an avalanche
[keep, IMA,mk,bk1,bk2,info,maxdifp,participationratio,Num_p, standardev] = discriminate(folder,En,ni,D,w,Cutoff,MinSep);

if (keep == 3)
    nfo = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d%05d.mat',folder,folder,En,En,ni);
    nfn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,ni);
    

        positions_exists = whos(matfile(nfo),'pxs','pys','Npf');
        
        if(length(positions_exists) == 3)
            load(nfo,'pxs','pys','Npf');
            track = 0;
            %save(nfn,'maxdifp','participationratio','Num_p', 'standardev','-append');
        else
            %% Particle centers
            pxs = zeros(3500,info.Numframe);
            pys = zeros(3500,info.Numframe);
            Npf = zeros(1,info.Numframe);
            surfacebestfitline = zeros(2,info.Numframe);
            %% Create background
            bg = bk1*0;
            bbg = bg;
            bbg2 = bbg;
            
            %% Create filter
            se = strel('disk',D+4);
            filtercof = 50;
            
            %% Parameters for Original image and ideal one
            
            % setup for ideal particle
            
            ss = 2*fix(D/2+4*w/2)-1;         % size of ideal particle image
            os = (ss-1)/2;                   % (size-1)/2 of ideal p
            [xx, yy] = ndgrid(-os:os,-os:os);  % ideal particle image grid
            rr = abs(xx+1i*yy);    % radial coordinate
            
        end
    
    
    %% Main loop over images
    for nframe = 1:info.Numframe
        
        if(track)
            
            im = IMA(:,:,nframe);
            bgmask = (imopen(im.*mk./bk1,se))>filtercof;
            normcoef = sum(sum(im.*bgmask))/sum(bgmask(:));
            
            im = im/normcoef;
            
            bg = max(bk1,im);
            bbg(im>0) = max(bk2(im>0),bg(im>0)-im(im>0));
            bbg2(bgmask) = bg(bgmask);
            bbg2(bgmask==0) = bbg(bgmask==0);
            
            
            low = 0.1;
            high = 0.95;
            sim = (clip((bg-im)./bbg2,low,high)-low)/(high-low);
            A = isnan(sim);
            sim(A) = 0;
            
            
            %Chi image ipf is the ideal image
            [ichi] = chiimg(sim,ipf(rr,D,w),[],[],'same');
            % find pixel accurate centers using chi-squared
            
            [~, py, px] = findpeaks(mk./ichi,mk,Cutoff,MinSep);  % find maxima
        else
            px = pxs(:,nframe);
            py = pys(:,nframe);
        end
        
        % Keep only insiders
        binsize = 15;
        [xs ix] = sort(ceil(px/binsize));
        bindex = [0 ;find(diff(xs)>0)];
        ys = py(ix);
        ya = zeros(1,length(bindex)-1);
        xa = ya;
        
        for bin = 1:length(bindex)-1;
            ya(bin) = min(ys(bindex(bin)+1:bindex(bin+1)));
            xa(bin) = binsize*(bin + xs(1)-1);
        end
        
        linearfit = polyfit(xa,ya,1);
        nvector = [-linearfit(1) 1]/(sqrt(linearfit(1)^2+1)); %normalvector to bestfit line
        lo = linearfit(2)*nvector(2); %Distance of bestfitline to origin.
        dpointtoline = px*nvector(1) + py*nvector(2) - lo;
        insiders = (dpointtoline > -D*7); %Points that are in the region of interest
        surfacebestfitline(:,nframe) = linearfit;
        
        %Get rid of outliers
        Npf(nframe) = sum(insiders);
        pxs(1:Npf(nframe),nframe) = px(insiders);
        pys(1:Npf(nframe),nframe) = py(insiders);
        
    end
    
    
    
    
    save(nfn,'git_version','pys','pxs','Npf','surfacebestfitline',...
        'maxdifp','participationratio','Num_p', 'standardev');
else
    save (sprintf('/aline%i/rotdrum%i/o%02d/no_avalanche_in_this_file%02d_%05d.mat',folder,folder,En,En,ni),...
        'git_version','maxdifp','participationratio','Num_p', 'standardev');
end





