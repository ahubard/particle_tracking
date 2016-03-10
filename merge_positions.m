function [x_ima, y_ima, x_non_ima, y_non_ima,overlap_size] = merge_positions(x_ima,y_ima,x_rot,y_rot)
% Find the overlaping regions of ima and rot

adjacentmatrix = adjacent(x_ima, y_ima, x_rot, y_rot, 7);
overlap_r = find(sum(adjacentmatrix));
i_non_overlap_r = setdiff((1:length(x_rot)),overlap_r);
overlap_i = find(sum(adjacentmatrix,2));
i_non_overlap_i = setdiff(1:length(x_ima),overlap_i);
overlap_size(1) = length(overlap_i);
overlap_size(2) = length(overlap_r);

% define overlaping positions
x_over_i = x_ima(overlap_i);
y_over_i = y_ima(overlap_i);
x_over_r = x_rot(overlap_r);
y_over_r =  y_rot(overlap_r);

% update image by non overlaping of original plus rotation
x_ima = [x_ima(i_non_overlap_i); x_rot(i_non_overlap_r)];
y_ima = [y_ima(i_non_overlap_i); y_rot(i_non_overlap_r)];

x_non_ima = x_rot(i_non_overlap_r);
y_non_ima = y_rot(i_non_overlap_r);

if(~isempty(overlap_i))
    
    % Assign closest trivial neighbours of overlaping region;
    [adjacentmatrix, tb1, tb2,~] = ...
        adjacent(x_over_i,y_over_i,x_over_r ,y_over_r,4.5);
    
    % Check overlaps
    overlap_r = find(sum(adjacentmatrix))';
    overlap_i = find(sum(adjacentmatrix,2));
    
    % Get the non trivial overlaps and keep one set;
    ntb1 = setdiff(overlap_i,tb1);
    ntb2 = setdiff(overlap_r,tb2);
    
    %**   %Get the non_overlaping
    %         x_ima = [x_ima; x_over_i(tb1)];
    %         y_ima = [y_ima; y_over_i(tb1)];
    
    x_ima = [x_ima; 1*x_over_i(tb1)+0*x_over_r(tb2)];
    y_ima = [y_ima; 1*y_over_i(tb1)+0*y_over_r(tb2)];
    
    i_non_r = setdiff((1:length(x_over_r)),overlap_r);
    i_non_i = setdiff(1:length(x_over_i),overlap_i);
    
    
    i1 = [ntb1; i_non_i'];
    i2 = [ntb2; i_non_r'];
    
    if(~isempty(i1) || ~isempty(i2))
        x1 = x_over_i(i1);
        y1 = y_over_i(i1);
        x2 = x_over_r(i2);
        y2 = y_over_r(i2);
        
        if(~isempty(i1))
            if(~isempty(i2))
                adjacentmatrix = adjacent(x1,y1,x2 ,y2,6);
                
                if(length(i1) >= length(i2))
                    non_ove = find(sum(adjacentmatrix) == 0);
                    x_extra = [x1; x2(non_ove)];
                    y_extra = [y1; y2(non_ove)];
                    x_non_ima = [x_non_ima; x2(non_ove)];
                    y_non_ima = [y_non_ima; y2(non_ove)];
                    
                else
                    non_ove = find(sum(adjacentmatrix,2) == 0);
                    x_extra = [x1(non_ove); x2];
                    y_extra = [y1(non_ove); y2];
                    x_non_ima = [x_non_ima; x2];
                    y_non_ima = [y_non_ima; y2];
                end
            else
                x_extra = x1;
                y_extra = y1;
            end
        else
            x_extra = x2;
            y_extra = y2;
            x_non_ima = [x_non_ima; x2];
            y_non_ima = [y_non_ima; y2];
        end
        x_ima = [x_ima; x_extra];
        y_ima = [y_ima; y_extra];
    end
end


%             plot(x_ima,y_ima,'.');axis('equal');
%             drawnow;
%             pause()
%**



%         if(~isempty(ntb1) || ~isempty(ntb2))
%             if(length(ntb1) < length(ntb2))
%                 x_extra = x_over_r(ntb2);
%                 y_extra = y_over_r(ntb2);
%             else
%                 x_extra = x_over_i(ntb1);
%                 y_extra = y_over_i(ntb1);
%             end
%         else
%             x_extra = [];
%             y_extra = [];
%         end

%Update list with weighted average of closest neighbours

%         x_ima = [x_ima; 0.8*x_over_i(tb1)+0.2*x_over_r(tb2); x_extra];
%         y_ima = [y_ima; 0.8*y_over_i(tb1)+0.2*y_over_r(tb2); y_extra];

%         %Get the non_overlaping
%         i_non_r = setdiff((1:length(x_over_r)),overlap_r);
%         i_non_i = setdiff(1:length(x_over_i),overlap_i);

%         if(~isempty(i_non_i) || ~isempty(i_non_r))
%             x_i = x_over_i(i_non_i);
%             y_i = y_over_i(i_non_i);
%             x_r = x_over_r(i_non_r);
%             y_r = y_over_r(i_non_r);
%             adjacentmatrix = adjacent(x_i,y_i,x_r ,y_r,6);
%             non_ove = find(sum(adjacentmatrix) == 0);
%             x_extra = [x_i; x_r(non_ove)];
%             y_extra = [y_i; y_r(non_ove)];
%          else
%             x_extra = [];
%             y_extra = [];
%         end
%
%             %Update list with weighted average of closest neighbours
%         x_ima = [x_ima; x_extra];
%         y_ima = [y_ima; y_extra];