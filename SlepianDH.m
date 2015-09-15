function [D]= SlepianDH(L_max,tt,pp,R_mask)
%SlepianDH generate the Hermitian Slepian D matrix for mesh region
% \int_R I_R(\unit{x}) Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})
% L_max is the highest non-zero degree
% tt and pp define the integration region R
% R_mask irregular region mask for I_R(\unit{x}) with 0-1 matrix on tt-pp grid 

if nargin < 4
	R_mask=0; % no mask; trapSphereMaskedR will bypass applying mask
end

N_tot=(L_max+1)^2;
D=zeros(N_tot,N_tot); % pre-allocate D
total=N_tot*(N_tot+1)/2; % number of elements to compute for Hermitian
count=0; progressbar('overall','inner loop')
for r=1:N_tot % rows of D
	if count~=0; progressbar([],0); end;
	l=floor(sqrt(r-1)); m=r-1-l*(l+1); % (7.40) here n=r-1
	f=sphereHarm(l,m,tt,pp);
	D(r,r)=trapSphereMaskedR(f.*conj(f),tt,pp,R_mask); % diagonal element
	count=count+1; progressbar(count/total,[])
	for c=r+1:N_tot % columns of D
		p=floor(sqrt(c-1)); q=c-1-l*(l+1); % (7.40) here n=c-1
		g=sphereHarm(p,q,tt,pp);
		D(r,c)=trapSphereMaskedR(g.*conj(f),tt,pp,R_mask); % (8.27)
		D(c,r)=conj(D(r,c)); % Hermitian symmetric element
		count=count+1; progressbar(count/total,(c-r)/(N_tot-r))
	end
	progressbar(count/total,[])
end
