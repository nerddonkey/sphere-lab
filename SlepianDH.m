function [D]= SlepianDH(L_max,tt,pp)
%SlepianDH generate the Hermitian Slepian D matrix for mesh region
% \int_R Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})
% L_max is the highest non-zero degree
% tt and pp define the integration region R
	N_tot=(L_max+1)^2;
	D=zeros(N_tot,N_tot); % pre-allocate D
	total=N_tot*(N_tot+1)/2; % number of elements to compute for Hermitian
	count=0; progressbar('overall','inner loop')
	for r=1:N_tot % rows
		progressbar([],0)
		l=floor(sqrt(r-1)); m=r-1-l*(l+1); % (7.40) n=r-1
		f=sphereHarm(l,m,tt,pp);
		D(r,r)=trapSphereR(f.*conj(f),tt,pp); % diagonal element
		count=count+1; progressbar(count/total,[])
		for c=r+1:N_tot % columns
			l=floor(sqrt(c-1)); m=c-1-l*(l+1); % (7.40) n=c-1
			g=sphereHarm(l,m,tt,pp);
			D(r,c)=trapSphereR(g.*conj(f),tt,pp);
			D(c,r)=conj(D(r,c)); % Hermitian symmetric element
			count=count+1; progressbar(count/total,(c-r)/(N_tot-r))
		end
		progressbar(count/total,[])
	end
end

% sphereHarm:
% Y_l^m = (-1)^m Q_l P_l^m(sch) e^m\phi, for m\geq0
% Y_l^{-m} = (-1)^m conj(Y_l^m), for m>0
% trapSphereR: integration of region R defined by tt and pp
% \int_R Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})
