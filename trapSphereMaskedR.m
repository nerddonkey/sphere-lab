function [I]=trapSphereMaskedR(f,tt,pp,R_mask)
%trapSphereMaskedR do \int_R m(x) f(x) ds(x) using 2D trapezoid rule integration
% tt and pp should be the separable mesh that covers region R
% an irregular region R is specified by the mask

if nargin<4
	R_mask=0;
end

fs=bsxfun(@times,f,sin(tt(:))); % apply sphere measure
if R_mask~=0
	fs=fs.*R_mask;
end
I=trapz(pp,trapz(tt,fs,1));
