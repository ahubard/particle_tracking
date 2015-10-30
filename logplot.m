
variable_number_id = 1;  %Gives which plot to create and their values



figure(1)
if( variable_number_id == 1)
    variable = -Size_Potential(ii_non_spaning);
    %variable = -Size_Potential(ii_whole);
    exponent = -1.17;
    normalization_coeff = .26;
    minnorm = 1;
elseif ( variable_number_id == 2)
    variable = T(ii_bound_left);
    exponent = -1.3;
    normalization_coeff = .12;
    minnorm = 1/60;
end

N_bin =31;
 xmin = log10(max(minnorm,min(variable)));
 xmax = log10(max(variable)+10*minnorm);
 x = logspace(xmin,xmax,N_bin);
%x = 1:3:124;
h = histogram(variable,x,'Normalization','pdf');
y = h.Values;
figure(2)
loglog(x(2:end),y1,'.-',x(2:end),y,'.-',x(2:end),normalization_coeff*x(2:end).^exponent);
loglog(x(2:end),y,'.-m');
    
