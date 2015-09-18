function [Ylm,theta,phi]=sphereHarm(l,m,tt,pp)
%sphereHarm spherical harmonic of degree l and order m
% see notes.tex/pdf for maths; seems robust for l to 2000

Qlm=(-1)^m*sqrt((2*l+1)/(8*pi)); % (-1)^m Q_l in note.tex/pdf

if isvector(tt) && isvector(pp)
	[theta,phi]=ndgrid(tt,pp); % create output mesh
	Sl=legendre(l,cos(tt),'sch');
	% Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 adjustment, see note.tex/pdf
	Slm=Sl(abs(m)+1,:)'; % pull out S_l^|m| row and transpose to column
	Ylm=Qlm*kron(ones(1,length(pp)),Slm).*exp(1j*abs(m)*phi);
elseif ~isvector(tt) && ~isvector(pp) && isequal(size(tt),size(pp))
	theta=tt; phi=pp; % pass input mesh to output
	S_max=numel(theta);
	for s=1:S_max % walk thru mesh as vector
		Sl=legendre(l,cos(theta(s)),'sch'); % vector in m, scalar in theta
		Slm=Sl(abs(m)+1,1); % pull out scalar S_l^|m| (scalar theta)
		Ylm(s)=Qlm*Slm.*exp(1j*abs(m)*phi(s)); % scalar
	end
else
	error('@@ sphereHarm: invalid inputs')
end

if m==0
	Ylm=sqrt(2)*Ylm; % m=0 adjustment, see note.tex/pdf
elseif m<0
	Ylm=(-1)^m*conj(Ylm); % m<0 adjustment, see note.tex/pdf
end
