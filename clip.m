function out=clip(in,lo,hi,ulo,uhi)
% clip Clips the values of in between lo and hi.
% Usage:: out=clip(in,lo,hi)
%
% Clips the values of in between lo and hi.

% revision history:
% 05/02/95 Mark D. Shattuck <mds>   clip.m

if(~exist('ulo','var')||isempty(ulo));
	ulo=lo;
end
if(~exist('uhi','var')||isempty(uhi));
	uhi=hi;
end

out=in;
out((out<lo))=ulo;
out((out>hi))=uhi;
