function  h = log2dhist(variable1,variable2,bins1,bins2)

%Creates a 2 dimensional histogram with the specified binning
h = zeros(length(bins1)-1,length(bins2)-1);
bin_area = diff(reshape(bins1,length(bins1),1))*diff(reshape(bins2,1,length(bins2)));

[x, ii] = sort(variable1);
y = variable2(ii);

lowx = find(x >= bins1(1),1,'first');
for ix = 2:length(bins1)
    hix = lowx-1+find(x(lowx:end)< bins1(ix),1,'last');
    if(isempty(hix))

         h(ix-1,:) = 0;
      
    else
        [yx, ~] = sort(y(lowx:hix));
        lowy = find(yx >= bins2(1),1,'first');
        for iy = 2:length(bins2)
            hiy = lowy-1+find(yx(lowy:end) < bins2(iy),1,'last');
            if(isempty(hiy))

                h(ix-1,iy-1) = 0;
            
            else
     
                h(ix-1,iy-1) = length(lowy:hiy);
                
                
                lowy = hiy+1;
                
            end 
        end
        lowx = hix+1;
    end
end

h = h./bin_area;

        
        
        
    
    
    
    
    

