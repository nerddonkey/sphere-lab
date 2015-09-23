function [Ylm,theta,phi]=sphHarmGrid(l,m,tv,pv)
%sphHarmGrid
%	computes the spherical harmonic of degree l and order m
%	on a theta-phi cartesian product sampling grid via
%	Schmidt semi-normalized Associated Legendre functions
%
%	integers l, m are the degree and order of the spherical harmonic
%	vectors tv and pv in radians define the locations of the samplig grid
%
%	see notes.tex/pdf for spherical harmonic maths; seems robust for l to 2000
%

[theta,phi]=ndgrid(tv,pv); % create output mesh
Qlm=(-1)^m*sqrt((2*l+1)/(8*pi)); % (-1)^m Q_l in note.tex/pdf
Sl=legendre(l,cos(tv),'sch');
Slm=Sl(abs(m)+1,:)'; % pull out S_l^|m| row and transpose to column
Ylm=Qlm*kron(ones(1,length(pv)),Slm).*exp(1j*m*phi);
if m==0
	Ylm=sqrt(2)*Ylm; % m=0 adjustment, see note.tex/pdf
elseif m<0
	Ylm=(-1)^m*Ylm; % m<0 adjustment, see note.tex/pdf
end
