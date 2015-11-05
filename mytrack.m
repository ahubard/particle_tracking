%Find the trayectories of the particles between files initial and final 
%using Hungaring  optimization algorithm assignmentoptimal


function  nfiles = mytrack(folder,En,NEn,initial,final,D,mk)


[git_version, ~] = evalc('system(''git describe --dirty --alway'')');
[px,py,NPF,initial] = stickfiles(folder,En,initial,final,D);
 maskfile = sprintf('/aline%i/rotdrum%i/o%i/mask%i.mat',folder,folder,En,En);
 if(~exist(maskfile,'file'))
     maskfile = sprintf('/aline%i/rotdrum%i/o%i/mask%i.mat',1,1,103,103);
 end
 load(maskfile,'mk');
nfiles = final-initial+1;

if (nfiles)
    
    t0 = 1;
    NY = size(mk,1);
    maxy = NY-2*D;
    maxdisptrivial = D/2;  %smaller than particle diameter
    maxdispnontrivial = 3*D/2;
    NPmax = max(NPF);
    NT = size(px,2);
    PX = zeros(NPmax,NT);
    PY = zeros(NPmax,NT);
    idtold = 1:NPmax;
    x1 = px(1:NPF(t0),t0);
    y1 = py(1:NPF(t0),t0) ;
    
    
    %% Loop over time
    for t1 = t0:NT-1
        
        t2=t1+1;
        
        
        %% Positions to track from frame 0 to frame 1
        x2 = px(1:NPF(t2),t2);
        y2 = py(1:NPF(t2),t2);
        NPT1 = length(x1);
        
        %Get adjacent matrix for distances smaller than maxdisp by using function
        %that bins the positions
        [~, trivialbondt1,trivialbondt2] = adjacent(x1,y1,x2,y2,maxdisptrivial);
        
        %Find particles that only have one bond from particle i(t1) to j(t2) so
        %conection is trivial.
        
        idt1 = trivialbondt1;
        idt2 = trivialbondt2;
        numtrackedt1 = length(idt1);
        numtrackedt2 = length(idt2);
        
        
        %Find the ids of the nontrivial bonds.
        if(length(trivialbondt1)<NPT1)
            
            mask = ones(1,NPT1);
            mask(trivialbondt1)=0;
            nontrivialt1 = find(mask);
            
            mask = ones(1,NPF(t2));
            mask(trivialbondt2)=0;
            nontrivialt2 = find(mask);
            
            x1n = x1(nontrivialt1);
            x2n = x2(nontrivialt2);
            y1n = y1(nontrivialt1);
            y2n = y2(nontrivialt2);
            
            %Check if particles are in the bottom and dont care about them
            nontrivialt1 = nontrivialt1(and(y1n<maxy,mk(sub2ind(size(mk),y1n,x1n))));
            nontrivialt2 = nontrivialt2(and(y2n<maxy,mk(sub2ind(size(mk),y2n,x2n))));
            x1n = x1(nontrivialt1);
            y1n = y1(nontrivialt1);
            x2n = x2(nontrivialt2);
            y2n = y2(nontrivialt2);
            
            if(length(x1n)<length(x2n))
                %warning('particle appeared')
                save(sprintf('/aline%i/rotdrum%i/o%02d/Warning1_%i_%i',folder,folder,En,NEn,t1),'x1n','y1n','x2n','y2n');
            end
            
            if(length(x1n)>length(x2n))
                %warning('particle disappeared')
                save(sprintf('/aline%i/rotdrum%i/o%02d/Warning2_%i_%i',folder,folder,En,NEn,t1),'x1n','y1n','x2n','y2n','PX','PY','idtold');
            end
            
            
            if(x1n) %Find the  posible bonds for particles that didnt have trivial bonds by allowing larger distances between them.
                %[length(x1n) length(x2n)]
                
                [adjacentmatrix, trivialbondt1, trivialbondt2, distMatrix] = adjacent(x1n,y1n,x2n,y2n,maxdispnontrivial);
                if (adjacentmatrix)
                    idt1(numtrackedt1+(1:length(trivialbondt1))) = nontrivialt1(trivialbondt1);
                    idt2(numtrackedt2+(1:length(trivialbondt2))) = nontrivialt2(trivialbondt2);
                    
                    
                    numtrackedt1 = length(idt1);
                    numtrackedt2 = length(idt2);
                    mask1 = ones(1,length(x1n));
                    mask1(trivialbondt1)=0;
                    nontrivialt1 = nontrivialt1((mask1>0));
                    nbnewtrackt1 = sum(mask1);
                    
                    mask2 = ones(1,length(x2n));
                    mask2(trivialbondt2)=0;
                    nontrivialt2 = nontrivialt2((mask2>0));
                    nbnewtrackt2 = sum(mask2);
                    
                    distMatrix = distMatrix(mask1>0,mask2>0);
                    distMatrix(adjacentmatrix(mask1>0,mask2>0)==0) = inf; %Set forbiden conections to infinity.
                    
                    [assignment, ~] = assignmentoptimal(distMatrix); %Call hungarian algorithm to find perfect matching
                    
                    if(nbnewtrackt1~=nbnewtrackt2)
                        
                        save(sprintf('/aline%i/rotdrum%i/o%02d/Warning3_%i_%i.mat',folder,folder,En,NEn,t1),'x1n','x2n','y1n','y2n','nontrivialt1','nontrivialt2');
                        if(nbnewtrackt1<nbnewtrackt2) %particles appeared that have no track
                            %warning('appeared particle stil here')
                            nbnewtrackt2 = nbnewtrackt1;
                        else %particle disappeared so assignment has zeros
                            %warning('its gone');
                            nontrivialt1 = nontrivialt1(assignment>0);
                            assignment = assignment(assignment>0);
                            nbnewtrackt1 =nbnewtrackt2;
                        end
                        
                        
                        
                        
                        
                        
                        
                    end
                    
                    if (sum(assignment == 0) > 0)
                        
                        save(sprintf('/aline%i/rotdrum%i/o%02d/Warning3_%i_%i.mat',folder,folder,En,NEn,t1),'x1n','x2n','y1n','y2n','nontrivialt1','nontrivialt2');
                        nontrivialt1 = nontrivialt1(assignment>0);
                        assignment = assignment(assignment>0);
                        nbnewtrackt1 =nbnewtrackt2;
                    end
                    
                    idt1(numtrackedt1+(1:nbnewtrackt1)) = nontrivialt1;
                    idt2(numtrackedt2+(1:nbnewtrackt2)) = nontrivialt2(assignment);
                end
            end
        end
        
        [idt1, orderidt1] = sort(idt1);
        idt2 = idt2(orderidt1);
        
        
        
        idt1 = idtold(idt1);
        % % %Save tracked particles
        %NPT2 = length(idt2);
        PX(idt1,t2) = x2(idt2);
        PY(idt1,t2) = y2(idt2);
        
        
        
        
        idtold = idt1;
        
        if( t1==1)
            PX((idt1),t1) = x1(idt1);
            PY((idt1),t1) = y1(idt1);
        end
        
        VX = PX(:,t2) - PX(:,t1);
        VY = PY(:,t2) - PY(:,t1);
        
        x1 = PX(idt1,t2) + VX(idt1) ;
        y1 = min(PY(idt1,t2) + VY(idt1),NY) ; % In case velocity make particle out of the boundary. 
        y1 = max(y1,1);  %in case y1 is negative. 
    end
    
    
    %% Save
    save(sprintf('/aline%i/rotdrum%i/o%02d/Tracked_%i.mat',folder,folder,En,NEn),'PX','PY','git_version');
    
else
    save(sprintf('/aline%i/rotdrum%i/o%02d/NoAvalanche_%i.mat',folder,folder,En,NEn),'nfiles','git_version');
end














