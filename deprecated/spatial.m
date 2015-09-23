function [F,theta,phi]=spatial(w,tt,pp,useProgressBar)
%SPATIAL performs the inverse spherical harmonic transform. The spherical harmonic
% coefficients are provided in vector w which is organized in the (l,m): (0,0)
% (1,-1) (1,0) (1,1) (2,-2) ... ordering.  Inputs tt and pp specify the theta and
% phi values over which the spatial function f is computed.  The coefficient vector
% (either row or column) will be zero padded to be a squared dimension long.
%
% The first case tt and pp are both vectors (either row or column) which results in
% the computation of a separable mesh, or rectangular grid, using built-in ndgrid.
% Function f, being the complex-valued weighted sum of complex spherical harmonics,
% is computed over this theta-phi mesh.  As well as f, the created mesh theta and
% phi is returned. Computations are efficient and vectorized.
%
% The second case tt and pp are both matrices of identical size representing a
% general, not necessarily rectangular, grid.  Theses are taken as the mesh and
% copied to the output mesh theta and phi.  Function f is computed individually at
% each point on the mesh.  Computations are essentially all scalar because we can't
% assume separability.
%
% Y_l^m = (-1)^m Q_l P_l^m(sch) e^m\phi, for m\geq0
% Y_l^{-m} = (-1)^m conj(Y_l^m), for m>0
%

if nargin < 4
	useProgressBar = 0;
end

% pad w with zeros if its length is not a square
L_max=ceil(sqrt(length(w)))-1;
N_tot=(L_max+1)^2;
if length(w) < N_tot % pad with zeros to make it a squared length
	w(N_tot)=0;
end

% SEPARABLE GRID CASE
% We exploit separability and the vector operations of matlab.
% The associated legendre function is done as a bank where rows
% are the order m and the columns are the vector of theta; this is
% for each degree l.  Similarly the complex exponentials are
% computed on a matrix mesh.  So the computation is pretty quick.
if isvector(tt) && isvector(pp)
	tt=tt(:)'; pp=pp(:)'; % ensure they are row vectors
	[theta,phi]=ndgrid(tt,pp); % define the mesh

	F=zeros(size(theta)); % mesh answer to be accumulated

	if useProgressBar~=0
		progressbar('spatial')
	end;

	for l=0:L_max
		Sl=legendre(l,cos(tt),'sch');
		Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 component, see note.tex/pdf
		Ql=sqrt((2*l+1)/(8*pi));
		for m=0:l
			n=l*(l+1)+m; % (7.39)
			n1=l*(l+1)-m; % (7.39)
			wlm=w(n+1); % the weight for degree l and order m
			wlm1=w(n1+1); % the weight for degree l and order m
			if wlm==0 && wlm1==0
				continue
         end
			Slm=Sl(m+1,:)'; % pull out S_l^m (evaluated on theta vector tt)
			Ylm=(-1)^m*Ql*kron(ones(size(pp)),Slm).*exp(1j*m*phi);
			F=F+wlm*Ylm; % m contribution
			if m~=0 % add the -m contribution
				F=F+wlm1*(-1)^m*conj(Ylm);
			end
		end
		if useProgressBar~=0
			progressbar((l+1)^2/(L_max+1)^2)
		end;
	end
	return;
end % SEPARABLE CASE

% GENERALLY NON-SEPARABLE GRID CASE
% The mesh is a rectangular configuration of nodes and in
% this case the rows and columns can be warped and curved
% in the theta phi plane.  Because the mesh is non-separable
% we need to compute for each point; so this can be slow.
if ~isvector(tt) && ~isvector(pp) && isequal(size(tt),size(pp))
	theta=tt; phi=pp; % pass input mesh to output
	S_max=numel(theta);

	F=zeros(size(theta)); % mesh answer to be accumulated

	if useProgressBar~=0
		progressbar('spatial non-sep')
	end;

	for s=1:S_max % walk thru matrix as vector
		for l=0:L_max
			Sl=legendre(l,cos(theta(s)),'sch');
			Sl(1,:)=Sl(1,:)*sqrt(2); % m=0 component
			Ql=sqrt((2*l+1)/(8*pi)); % Y_l^m = (-1)^m NSl P_l^m(sch) e^m\phi
			for m=0:l
				n=l*(l+1)+m; % (7.39)
				n1=l*(l+1)-m; % (7.39)
				wlm=w(n+1); % the weight for degree l and order m
				wlm1=w(n1+1); % the weight for degree l and order m
				if wlm==0 && wlm1==0
					continue
				end
				Slm=Sl(m+1,1); % pull out scalar S_l^m (evaluated at single theta)
				Ylm=(-1)^m*Ql*Slm.*exp(1j*m*phi(s)); % scalar
				F(s)=F(s)+wlm*Ylm;
				if m~=0
					F(s)=F(s)+wlm1*(-1)^m*conj(Ylm);
				end
			end
		end
		if useProgressBar~=0
			progressbar(s/S_max);
		end;
	end
	return;
end % GENERALLY NON-SEPARABLE GRID CASE

error('weird mesh specification')
