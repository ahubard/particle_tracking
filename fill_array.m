function  C = fill_array(nv,var1,var2, var3,var4,var5,var6,var7)

%Fills an array with vectors of different size filling the rest with zeros.
%There is probably an easier and smarter way. 

l = zeros(1,nv);

%Find sizes of variables
l(1) = length(var1);
l(2) = length(var2);
l(3) = length(var3);
l(4) = length(var4);
l(5) = length(var5);
l(6) = length(var6);
l(7) = length(var7);

array_length = max(l);
length_diff = array_length-l;
%fill cells with max seven arrays

C = zeros(nv,array_length);
C(1,1:l(1)) = var1(1:l(1));
C(2,1:l(2)) = var2(1:l(2));
C(3,1:l(3)) = var3(1:l(3));
C(4,1:l(4)) = var4(1:l(4));
C(5,1:l(5)) = var5(1:l(5));
C(6,1:l(6)) = var6(1:l(6));
C(7,1:l(7)) = var7(1:l(7));






