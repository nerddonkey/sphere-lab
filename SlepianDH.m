function [D]= SlepianDH(L_max,tt,pp)
%SlepianDH generate the Hermitian Slepian D matrix for mesh region
	N_max=(L_max+1)^2;
	D=zeros(N_max,N_max); % pre-allocate D
	total=N_max*(N_max+1)/2; % number of elements to compute for Hermitian
	count=0; progressbar('overall','inner loop')
	for r=1:N_max % rows
		progressbar([],0)
		if 0 % runs ISHT with delta function for weights
			f=spatial_shn(r,tt,pp); % Y_r or Y_l^m in (8.27)
		else
			l=floor(sqrt(r-1)); m=r-1-l*(l+1); % (7.40)
			f=sphereHarm(l,m,tt,pp);
		end
		D(r,r)=trapSphereR(f.*conj(f),tt,pp); % diagonal element
		count=count+1; progressbar(count/total,[])
		for c=r+1:N_max % columns
			if 0
				g=spatial_shn(c,tt,pp); % Y_c or Y_p^q in (8.27)
			else
				l=floor(sqrt(c-1)); m=c-1-l*(l+1); % (7.40)
				g=sphereHarm(l,m,tt,pp);
			end
			D(r,c)=trapSphereR(g.*conj(f),tt,pp);
			D(c,r)=conj(D(r,c)); % Hermitian symmetry
			count=count+1; progressbar(count/total,(c-r)/(N_max-r))
		end
		progressbar(count/total,[])
	end
end

% Y_l^m = (-1)^m Q_l P_l^m(sch) e^m\phi, for m\geq0
% Y_l^{-m} = (-1)^m conj(Y_l^m), for m>0
