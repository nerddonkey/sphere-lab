function [I]=trapSphereR(f,tt,pp)
%trapSphereR do \int_R f(x) ds(x) using 2D trapezoid rule integration
	fs=bsxfun(@times,f,sin(tt(:)));
	I=trapz(pp,trapz(tt,fs,1));
end
