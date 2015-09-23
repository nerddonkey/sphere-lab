function [YB,theta,phi]=sphHarmBankGrid(l,tv,pv)
%sphHarmBankGrid
%	computes the spherical harmonic of degree l and orders m=0:l
%	on a theta-phi cartesian product sampling grid via
%	Schmidt semi-normalized Associated Legendre functions
%
%	integer l is the degree of the spherical harmonic
%	vectors tv and pv in radians define the locations of the samplig grid
%
%	see notes.tex/pdf for spherical harmonic maths; seems robust for l to 2000
%

[theta,phi]=ndgrid(tv,pv); % create output mesh
Sl=legendre(l,cos(tv),'sch');
Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
Ql=sqrt((2*l+1)/(8*pi)); % Q_l in note.tex/pdf
for m=0:l
	Slm=Sl(m+1,:)'; % pull out S_l^m row and transpose to column
	YB{m+1}=(-1)^m*Ql*kron(ones(1,length(pv)),Slm).*exp(1j*m*phi); % cell array
end


function [Ylm,theta,phi]=sphHarmBankmGrid(l,m,tv,pv)
%sphHarmBankmGrid spherical harmonic of degree l and order m

[YB,theta,phi]=sphHarmBankGrid(l,tv,pv);
Ylm=YB{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end
