function [f,theta,phi]=spatial(w,tt,pp)
% SPATIAL performs the inverse spherical harmonic transform. The spherical harmonic
% coefficients are provided in vector w which is organized in the (l,m): (0,0)
% (1,-1) (1,0) (1,1) (2,-2) ... ordering.  Inputs tt and pp specify the theta and
% phi values over which the spatial function f is computed.  The coefficient vector
% (either row or column) should be L^2 long, otherwise it is implicitly truncated
% back to a shorter vector of square dimension.
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
% test for errors
	if ~isvector(w)
	  error('weight is not a vector')
	end

	L_max=floor(sqrt(numel(w)))-1;

	% SEPARABLE GRID CASE
	% We exploit separability and the vector operations of matlab.
	% The associated legendre function is done as a bank where rows
	% are the order m and the columns are the vector of theta; this is
	% for each degree l.  Similarly the complex exponentials are
	% computed on a matrix mesh.  So the computation is pretty quick.
	if isvector(tt) && isvector(pp)
		tt=tt(:)'; pp=pp(:)'; % ensure they are row vectors
		[theta,phi]=ndgrid(tt,pp); % define the mesh

		f=zeros(size(theta)); % mesh answer to be accumulated

		disp(sprintf('\n* Progress\n'))
		tic

		for l=[0:L_max]
			Sl=legendre(l,cos(tt),'sch'); % [S_l^0; ...; S_l^l](l,tt)
			NSl=sqrt((2*l+1)/(8*pi)); % Y_l^m = (-1)^m NSl P_l^m(sch) e^m\phi
			for m=[0:l]
				n=l*(l+1)+m; % (7.39)
				n1=l*(l+1)-m; % (7.39)
				wlm=w(n+1); % the weight for degree l and order m
				wlm1=w(n1+1); % the weight for degree l and order m
				if wlm==0 && wlm1==0; continue; end
				Slm=Sl(m+1,:)'; % pull out S_l^m (evaluated on theta vector tt)

				Ylm=(-1)^m*NSl*kron(ones(size(pp)),Slm).*exp(1j*m*phi);
				if m==0
					f=f+wlm*Ylm;
				else
					f=f+wlm*Ylm+wlm1*(-1)^m*conj(Ylm);
				end
			end
			disp(sprintf('%6d in seconds: %.2f', l, toc))
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

		f=zeros(size(theta)); % mesh answer to be accumulated

		for s=[1:numel(theta)] % walk thru matrix as vector
			for l=[0:L_max]
%	 			Pl=legendre3(l,cos(theta(s)));
%	 			Nl0=sqrt((2*l+1)/(4*pi)); % (7.20) for m=0
				Sl=legendre(l,cos(theta(s)),'sch');
				NSl=sqrt((2*l+1)/(8*pi)); % Y_l^m = (-1)^m NSl P_l^m(sch) e^m\phi
				for m=[0:l]
%	 				n=l*(l+1)+m; % (7.39)
%	 				wlm=w(n+1); % the weight for degree l and order m
% 	 				if wlm==0
% 	 					continue; % skip to save unnecessary computation
% 	 				end
% 	 				Plm=Pl(l+m+1,1); % pull out scalar P_l^m(sTheta)
% 					Nlm=Nl0*sqrt(factorial(l-m)/factorial(l+m)); % normalization (7.20)
% 	 				Ylm=Nlm*Plm*exp(1j*m*phi(s)); % scalar
% 	 				f(s)=f(s)+wlm*Ylm; % accumulate mesh answer here per mesh point
					n=l*(l+1)+m; % (7.39)
					n1=l*(l+1)-m; % (7.39)
					wlm=w(n+1); % the weight for degree l and order m
					wlm1=w(n1+1); % the weight for degree l and order m
					if wlm==0 && wlm1==0; continue; end
					Slm=Sl(m+1,:)'; % pull out S_l^m (evaluated on theta vector tt)

					Ylm=(-1)^m*NSl*Slm.*exp(1j*m*phi(s)); % scalar
					if m==0
						f(s)=f(s)+wlm*Ylm;
					else
						f(s)=f(s)+wlm*Ylm+wlm1*(-1)^m*conj(Ylm);
					end
				end
			end
		end
		return;
	end % GENERALLY NON-SEPARABLE GRID CASE
	error('weird mesh specification')
end
