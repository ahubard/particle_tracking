function ipi=ipr(cr,D,d,w)
% ipf    Calculate ideal particle image for rings. 
% Usage: ipi=ipf(cr,D,w)
%
% Calculates an ideal particle image ipi.  The particle has big diameter D+d
% and small one D-d.
% width parameter w.  2w is the width of 76% of the fall off. For the dot


ipi=(ipf(cr,D+d,w)-ipf(cr,D-d,w))/(ipf(D/2,D+d,w)-ipf(D/2,D-d,w));

