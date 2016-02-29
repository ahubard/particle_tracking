function [xn, yn] = rot_me(theta,x,y)

%Rotate the 2 dimentional vector (x,y) theta degrees around vector(0,0), x
%and y must have the same size and be vectors or scalars
if(length(theta) ~= length(x))
    if(length(theta) == 1)
        theta = repmat(theta,1,length(x));
    else
        error(' size of theta must be one or the same as x');
    end
end


s = sin(theta(:));
c = cos(theta(:));

xn = x(:).*c-y(:).*s;
yn = x(:).*s+y(:).*c;