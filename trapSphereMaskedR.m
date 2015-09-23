function intR=trapSphereMaskedR(f,tv,pv,R_mask)
%trapSphereMaskedR do \int_R m(x) f(x) ds(x) using 2D trapezoid rule integration
%
%  tv and pv should be the separable mesh that covers region R
%  an irregular region R is specified by the mask R_mask
%

%% Default to no mask
if nargin<4
	R_mask=0;
end

%% Apply sphere measure
fs=bsxfun(@times,f,sin(tv(:)));

%% Apply the mask
if R_mask~=0
	fs=fs.*R_mask;
end

%% Numerically integrate
intR=trapz(pv,trapz(tv,fs,1));
