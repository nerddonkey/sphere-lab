function [Ylm,theta,phi]=sphereHarmBankm(l,m,tt,pp)
%sphereHarmBankm
%	helper function to extract a single spherical harmonic
%	from a call to sphereHarmBank
%
%	computes the spherical harmonic of degree l and order m
%	on a theta-phi cartesian product sampling grid via
%	Schmidt semi-normalized Associated Legendre functions
%
%	integers l, m are the degree and order of the spherical harmonic
%	vectors tv and pv in radians define the locations of the samplig grid
%
%	see notes.tex/pdf for spherical harmonic maths; seems robust for l to 2000
%

[YB,theta,phi]=sphereHarmBank(l,tt,pp);
Ylm=YB{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end
